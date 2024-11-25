#include "Poisson2D.hpp"

void
Poisson2D::setup()
{
  std::cout << "===============================================" << std::endl;

  // Create the mesh.
  {
    std::cout << "Initializing the mesh from " << mesh_file_name << std::endl;

    // Or read the mesh from file:
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);

    std::ifstream mesh_file(mesh_file_name);
    grid_in.read_msh(mesh_file);

    std::cout << "  Number of elements = " << mesh.n_active_cells()
              << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    std::cout << "Initializing the finite element space" << std::endl;

    // Construct the finite element object. Notice that we use the FE_SimplexP
    // class here, that is suitable for triangular (or tetrahedral) meshes.
    fe = std::make_unique<FE_SimplexP<dim>>(r);

    std::cout << "  Degree                     = " << fe->degree << std::endl;
    std::cout << "  DoFs per cell              = " << fe->dofs_per_cell
              << std::endl;

    // Construct the quadrature formula of the appopriate degree of exactness.
    // This formula integrates exactly the mass matrix terms (i.e. products of
    // basis functions).
    quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

    std::cout << "  Quadrature points per cell = " << quadrature->size()
              << std::endl;

    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(r + 1);

    std::cout << "  Quadrature points per boundary cell = "
              << quadrature_boundary->size() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    std::cout << "Initializing the DoF handler" << std::endl;

    // Initialize the DoF handler with the mesh we constructed.
    dof_handler.reinit(mesh);

    // "Distribute" the degrees of freedom. For a given finite element space,
    // initializes info on the control variables (how many they are, where
    // they are collocated, their "global indices", ...).
    dof_handler.distribute_dofs(*fe);

    std::cout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  std::cout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    std::cout << "Initializing the linear system" << std::endl;

    // We first initialize a "sparsity pattern", i.e. a data structure that
    // indicates which entries of the matrix are zero and which are different
    // from zero. To do so, we construct first a DynamicSparsityPattern (a
    // sparsity pattern stored in a memory- and access-inefficient way, but
    // fast to write) and then convert it to a SparsityPattern (which is more
    // efficient, but cannot be modified).
    std::cout << "  Initializing the sparsity pattern" << std::endl;
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // Then, we use the sparsity pattern to initialize the system matrix
    std::cout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity_pattern);

    // Finally, we initialize the right-hand side and solution vectors.
    std::cout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(dof_handler.n_dofs());
    std::cout << "  Initializing the solution vector" << std::endl;
    solution.reinit(dof_handler.n_dofs());
  }
}

void
Poisson2D::assemble()
{
  std::cout << "===============================================" << std::endl;

  std::cout << "  Assembling the linear system" << std::endl;

  // Number of local DoFs for each element.
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  // Number of quadrature points for each element.
  const unsigned int n_q = quadrature->size();

  // FEValues instance. This object allows to compute basis functions, their
  // derivatives, the reference-to-current element mapping and its
  // derivatives on all quadrature points of all elements.
  FEValues<dim> fe_values(
    *fe,
    *quadrature,
    // Here we specify what quantities we need FEValues to compute on
    // quadrature points. For our test, we need:
    // - the values of shape functions (update_values);
    // - the derivative of shape functions (update_gradients);
    // - the position of quadrature points (update_quadrature_points);
    // - the product J_c(x_q)*w_q (update_JxW_values).
    update_values | update_gradients | update_quadrature_points |
      update_JxW_values);

  // Since we need to compute integrals on the boundary for Neumann conditions,
  // we also need a FEValues object to compute quantities on boundary edges
  // (faces).
  FEFaceValues<dim> fe_values_boundary(*fe,
                                       *quadrature_boundary,
                                       update_values |
                                         update_quadrature_points |
                                         update_JxW_values);

  // Local matrix and right-hand side vector. We will overwrite them for
  // each element within the loop.
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // We will use this vector to store the global indices of the DoFs of the
  // current element within the loop.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // Reset the global matrix and vector, just in case.
  system_matrix = 0.0;
  system_rhs    = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // Reinitialize the FEValues object on current element. This
      // precomputes all the quantities we requested when constructing
      // FEValues (see the update_* flags above) for all quadrature nodes of
      // the current cell.
      fe_values.reinit(cell);

      // We reset the cell matrix and vector (discarding any leftovers from
      // previous element).
      cell_matrix = 0.0;
      cell_rhs    = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          // Here we assemble the local contribution for current cell and
          // current quadrature point, filling the local matrix and vector.

          // Here we iterate over *local* DoF indices.
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  // FEValues::shape_grad(i, q) returns the gradient of the i-th
                  // basis function at the q-th quadrature node, already mapped
                  // on the physical element: we don't have to deal with the
                  // mapping, it's all hidden inside FEValues.
                  cell_matrix(i, j) += diffusion_coefficient.value(
                                         fe_values.quadrature_point(q)) // mu(x)
                                       * fe_values.shape_grad(i, q)     // (I)
                                       * fe_values.shape_grad(j, q)     // (II)
                                       * fe_values.JxW(q);              // (III)
                }

              cell_rhs(i) += forcing_term.value(fe_values.quadrature_point(q)) *
                             fe_values.shape_value(i, q) * fe_values.JxW(q);
            }
        }

      // If the cell is adjacent to the boundary...
      if (cell->at_boundary())
        {
          // ...we loop over its edges (referred to as faces in the deal.II
          // jargon).
          for (unsigned int face_number = 0; face_number < cell->n_faces();
               ++face_number)
            {
              // If current face lies on the boundary, and its boundary ID (or
              // tag) is that of one of the Neumann boundaries, we assemble the
              // boundary integral.
              if (cell->face(face_number)->at_boundary() &&
                  (cell->face(face_number)->boundary_id() == 2 ||
                   cell->face(face_number)->boundary_id() == 3))
                {
                  fe_values_boundary.reinit(cell, face_number);

                  for (unsigned int q = 0; q < quadrature_boundary->size(); ++q)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      cell_rhs(i) +=
                        function_h.value(
                          fe_values_boundary.quadrature_point(q)) * // h(xq)
                        fe_values_boundary.shape_value(i, q) *      // v(xq)
                        fe_values_boundary.JxW(q);                  // Jq wq
                }
            }
        }


      // At this point the local matrix and vector are constructed: we
      // need to sum them into the global matrix and vector. To this end,
      // we need to retrieve the global indices of the DoFs of current
      // cell.
      cell->get_dof_indices(dof_indices);

      // Then, we add the local matrix and vector into the corresponding
      // positions of the global matrix and vector.
      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
    }

  // Boundary conditions.
  {
    // We construct a map that stores, for each DoF corresponding to a Dirichlet
    // condition, the corresponding value. E.g., if the Dirichlet condition is
    // u_i = b, the map will contain the pair (i, b).
    std::map<types::global_dof_index, double> boundary_values;

    // Then, we build a map that, for each boundary tag, stores the
    // corresponding boundary function.
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    boundary_functions[0] = &function_g;
    boundary_functions[1] = &function_g;

    // interpolate_boundary_values fills the boundary_values map.
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    // Finally, we modify the linear system to apply the boundary conditions.
    // This replaces the equations for the boundary DoFs with the corresponding
    // u_i = 0 equations.
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}

void
Poisson2D::solve()
{
  std::cout << "===============================================" << std::endl;

  // Here we specify the maximum number of iterations of the iterative solver,
  // and its tolerance.
  SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());

  // Since the system matrix is symmetric and positive definite, we solve the
  // system using the conjugate gradient method.
  SolverCG<Vector<double>> solver(solver_control);

  std::cout << "  Solving the linear system" << std::endl;
  // We don't use any preconditioner for now, so we pass the identity matrix as
  // preconditioner.
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  std::cout << "  " << solver_control.last_step() << " CG iterations"
            << std::endl;
}

void
Poisson2D::output() const
{
  std::cout << "===============================================" << std::endl;

  // The DataOut class manages writing the results to a file.
  DataOut<dim> data_out;

  // It can write multiple variables (defined on the same mesh) to a single
  // file. Each of them can be added by calling add_data_vector, passing the
  // associated DoFHandler and a name.
  data_out.add_data_vector(dof_handler, solution, "solution");

  // Once all vectors have been inserted, call build_patches to finalize the
  // DataOut object, preparing it for writing to file.
  data_out.build_patches();

  // Then, use one of the many write_* methods to write the file in an
  // appropriate format.
  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string           output_file_name =
    "output-" + mesh_path.stem().string() + ".vtk";
  std::ofstream output_file(output_file_name);
  data_out.write_vtk(output_file);

  std::cout << "Output written to " << output_file_name << std::endl;

  std::cout << "===============================================" << std::endl;
}