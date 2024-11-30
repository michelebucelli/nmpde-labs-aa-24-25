#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson2D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  const std::string mesh_file_name =
    "../../examples/gmsh/mesh-square-h0.100000.msh";
  const unsigned int r = 1;

  Poisson2D problem(mesh_file_name, r);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}