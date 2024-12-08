#include "NonLinearDiffusion.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name = "../mesh/mesh-cube-20.msh";
  const unsigned int degree         = 1;

  NonLinearDiffusion problem(mesh_file_name, degree);

  problem.setup();
  problem.solve_newton();
  problem.output();

  return 0;
}