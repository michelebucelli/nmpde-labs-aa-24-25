#include <iostream>

#include "Poisson1D.hpp"

// Main function.
int
main(int /*argc*/, char * /*argv*/[])
{
  const unsigned int N = 39;
  const unsigned int r = 1;

  Poisson1D problem(N, r);

  problem.setup();
  problem.assemble();
  problem.solve();
  problem.output();

  return 0;
}
