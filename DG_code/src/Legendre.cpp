#include "Legendre.hpp"

#include <iostream>

namespace PolyDG
{

Real legendre(unsigned n, Real x)
{
  switch(n)
  {
    case 0:
      return 1.0;

    case 1:
      return x;

    case 2:
      return 1.5 * pow(x, 2) - 0.5;

    case 3:
      return x * (2.5 * pow(x, 2) - 1.5);

    case 4:
      return (35.0 * pow(x, 4) - 30.0 * pow(x,2) + 3.0) / 8.0;

    case 5:
      return (63.0 * pow(x, 5) - 70.0 * pow(x, 3) + 15.0 * x) / 8.0;

    case 6:
      return (231.0 * pow(x, 6) - 315.0 * pow(x, 4) + 105.0 * pow(x, 2) - 5.0) / 16.0;

    default:
      std::cerr << "The required degree is not implemented." << std::endl;
      exit(1);
  }
}

Real legendreDer(unsigned n, Real x)
{
  switch(n)
  {
    case 0:
      return 0.0;

    case 1:
      return 1.0;

    case 2:
      return 3.0 * x;

    case 3:
      return 7.5 * pow(x, 2) - 1.5;

    case 4:
      return (140.0 * pow(x, 3) - 60.0 * x) / 8.0;

    case 5:
      return (315.0 * pow(x, 4) - 210.0 * pow(x, 2) + 15.0) / 8.0;

    case 6:
      return (1386.0 * pow(x, 5) - 1260.0 * pow(x, 3) + 210.0 * x) / 16.0;

    default:
      std::cerr << "The required degree is not implemented." << std::endl;
      exit(1);
  }
}

} // namespace PolyDG
