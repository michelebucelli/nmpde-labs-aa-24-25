#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson3D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  const std::string  mesh_filename = "../mesh/mesh-cube-10.msh";
  const unsigned int degree        = 1;

  Poisson3D problem(mesh_filename, degree);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}