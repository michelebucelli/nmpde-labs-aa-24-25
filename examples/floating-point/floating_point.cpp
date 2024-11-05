#include <iostream>
#include <limits>
#include <vector>

int
main(int /*argc*/, char ** /*argv*/)
{
  using NL = std::numeric_limits<double>;

  // Check whether our implementation is IEC60559-compliant.
  std::cout << "is_iec559 = " << NL::is_iec559 << std::endl;

  // Retrieve tha minimum and maximum values that can be represented.
  std::cout << "min       = " << NL::min() << std::endl;
  std::cout << "max       = " << NL::max() << std::endl;

  // Retrieve the machine epsilon.
  std::cout << "epsilon   = " << NL::epsilon() << std::endl;

  std::cout << "========================" << std::endl;

  // Compute a = 1 - 3 * (4 / 3 - 1): it is zero in exact arithmetic, but it is
  // not zero in floating point arithmetic.
  std::cout << "a         = " << (1.0 - 3.0 * (4.0 / 3.0 - 1.0)) << std::endl;
  std::cout << "========================" << std::endl;

  // Compute f(x) = ((1 + x) - 1) / x for several (small) values of x, and
  // compute the relative error against the exact value f(x) = 1.
  std::vector<double> x_vals = {1e-10, 1e-14, 1e-15, 1e-16};

  std::cout << " x\t    err(x) " << std::endl;
  std::cout << "------------------------" << std::endl;

  for (const auto &x : x_vals)
    {
      const double f_of_x = ((1.0 + x) - 1.0) / x;
      const double error  = std::abs(f_of_x - 1.0) * 100;

      std::cout << " " << x << "\t    " << error << std::endl;
    }

  return 0;
}