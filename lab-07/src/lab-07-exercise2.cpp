#include "LinearElasticity.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name = "../mesh/mesh-cube-10.msh";
  const unsigned int degree         = 1;

  LinearElasticity problem(mesh_file_name, degree);

  problem.setup();
  problem.assemble_system();
  problem.solve_system();
  problem.output();

  return 0;
}