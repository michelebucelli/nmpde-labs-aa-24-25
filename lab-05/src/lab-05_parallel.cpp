#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D_parallel.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  // This object calls MPI_Init when it is constructed, and MPI_Finalize when it
  // is destroyed. It also initializes several other libraries bundled with
  // dealii (e.g. p4est, PETSc, ...).
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_filename = "../mesh/mesh-cube-40.msh";
  const unsigned int degree        = 1;

  Poisson3DParallel problem(mesh_filename, degree);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}