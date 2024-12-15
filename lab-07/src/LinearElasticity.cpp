#include "LinearElasticity.hpp"

void
LinearElasticity::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // First we read the mesh from file into a serial (i.e. not parallel)
    // triangulation.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream grid_in_file(mesh_file_name);
      grid_in.read_msh(grid_in_file);
    }

    // Then, we copy the triangulation into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);
      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }

    // Notice that we write here the number of *global* active cells (across all
    // processes).
    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space. This is the same as in serial codes.
  {
    pcout << "Initializing the finite element space" << std::endl;

    // To construct a vector-valued finite element space, we use the FESystem
    // class. It is still derived from FiniteElement.
    FE_SimplexP<dim> fe_scalar(r);
    fe = std::make_unique<FESystem<dim>>(fe_scalar, dim);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We retrieve the set of locally owned DoFs, which will be useful when
    // initializing linear algebra classes.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // To initialize the sparsity pattern, we use Trilinos' class, that manages
    // some of the inter-process communication.
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);

    // After initialization, we need to call compress, so that all process
    // retrieve the information they need for the rows they own (i.e. the rows
    // corresponding to locally owned DoFs).
    sparsity.compress();

    // Then, we use the sparsity pattern to initialize the system matrix. Since
    // the sparsity pattern is partitioned by row, so will the matrix.
    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity);

    // Finally, we initialize the right-hand side and solution vectors.
    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
LinearElasticity::assemble_system()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;

  // This class allows us to access vector-valued shape functions, so that we
  // don't have to worry about dealing with their components, but we can
  // directly use the vectorial form of the weak formulation.
  FEValuesExtractors::Vector displacement(0);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          // Evaluate coefficients on this quadrature node.
          const double mu_loc     = mu.value(fe_values.quadrature_point(q));
          const double lambda_loc = lambda.value(fe_values.quadrature_point(q));

          // Evaluate the forcing term on this quadrature node.
          Vector<double> f_loc(dim);
          forcing_term.vector_value(fe_values.quadrature_point(q), f_loc);

          // Convert the forcing term to a tensor (so that we can use it with
          // scalar_product).
          Tensor<1, dim> f_loc_tensor;
          for (unsigned int i = 0; i < dim; ++i)
            f_loc_tensor[i] = f_loc[i];

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  cell_matrix(i, j) +=
                    (mu_loc *
                       scalar_product(fe_values[displacement].gradient(j, q),
                                      fe_values[displacement].gradient(i, q)) +
                     lambda_loc * fe_values[displacement].divergence(j, q) *
                       fe_values[displacement].divergence(i, q)) *
                    fe_values.JxW(q);
                }

              cell_rhs(i) +=
                scalar_product(f_loc_tensor,
                               fe_values[displacement].value(i, q)) *
                fe_values.JxW(q);
            }
        }

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // Boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    boundary_functions[0] = &function_g;
    boundary_functions[1] = &function_g;

    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution_owned, system_rhs, true);
  }
}

void
LinearElasticity::solve_system()
{
  pcout << "===============================================" << std::endl;
  pcout << "Solving the system" << std::endl;

  SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());

  SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR      preconditioner;
  preconditioner.initialize(
    system_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " CG iterations" << std::endl;

  solution = solution_owned;
}

void
LinearElasticity::output() const
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  // By passing these two additional arguments to add_data_vector, we specify
  // that the three components of the solution are actually the three components
  // of a vector, so that the visualization program can take that into account.
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  std::vector<std::string> solution_names(dim, "u");

  data_out.add_data_vector(dof_handler,
                           solution,
                           solution_names,
                           data_component_interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-linearelasticity";
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << "." << std::endl;

  pcout << "===============================================" << std::endl;
}