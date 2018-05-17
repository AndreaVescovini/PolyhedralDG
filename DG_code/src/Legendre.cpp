#include "Legendre.hpp"

#include <stdexcept>

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

    case 7:
      return (429.0 * pow(x, 7) - 693.0 * pow(x, 5) + 315.0 * pow(x, 3) - 35.0 * x) / 16.0;

    case 8:
      return (6435.0 * pow(x, 8) - 12012.0 * pow(x, 6) + 6930.0 * pow(x, 4) - 1260.0 * pow(x, 2) + 35.0) / 128.0;

    default:
      throw std::domain_error("The required degree for the FeSpace is not implemented.");
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

    case 7:
      return (3003.0 * pow(x, 6) - 3465.0 * pow(x, 4) + 945.0 * pow(x, 2) - 35.0) / 16.0;

    case 8:
      return (51480.0 * pow(x, 7) - 72072.0 * pow(x, 5) + 27720.0 * pow(x, 3) - 2520.0 * pow(x, 1))/ 128.0;

    default:
      throw std::domain_error("The required degree for the FeSpace is not implemented.");
  }
}

} // namespace PolyDG
