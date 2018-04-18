#include "Legendre.hpp"

// #include <boost/math/special_functions/legendre.hpp>

#include <iostream>

namespace basis
{

std::array<PolyDG::Real, 2> legendre(unsigned n, PolyDG::Real x)
{
  switch(n)
  {
    case 0:
      return {1., 0};

    case 1:
      return {x, 1.};

    case 2:
      return {1.5 * pow(x, 2) - 0.5,
              3. * x};

    case 3:
      return {x * (2.5 * pow(x, 2) - 1.5),
              7.5 * pow(x, 2) - 3.};

    case 4:
      return {(35. * pow(x, 4) - 30. * pow(x,2) + 3.) / 8.,
              (140. * pow(x, 3) - 60. * x) / 8.};

    case 5:
      return {(63. * pow(x, 5) - 70. * pow(x, 3) + 15. * x) / 8.,
              (315. * pow(x, 4) - 210. * pow(x, 2) + 15.) / 8.};

    case 6:
      return {(231. * pow(x, 6) - 315. * pow(x, 4) + 105. * pow(x, 2) - 5.) / 16.,
              (1386. * pow(x, 5) - 1260. * pow(x, 3) + 210. * x) / 16.};

    default:
    // return {boost::math::legendre_p(n, x), boost::math::legendre_p_prime(n, x)};
      std::cerr << "Boost is needed." << std::endl;
      return {0, 0};

  }

}

} // namespace basis
