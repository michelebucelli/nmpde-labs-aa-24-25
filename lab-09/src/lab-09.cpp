#include "Stokes.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name  = "../mesh/mesh-step-5.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  Stokes problem(mesh_file_name, degree_velocity, degree_pressure);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}