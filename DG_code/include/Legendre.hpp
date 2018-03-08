#ifndef _LEGENDRE_HPP_
#define _LEGENDRE_HPP_

#include "geom.hpp"
#include <array>

namespace basis
{

template<class T>
constexpr T pow(const T& x, unsigned n);

// Function that return the value of the legendre polynomial of order n at the
// point x belonging to [-1, 1] and the value of the derivative at the same point.
std::array<geom::real, 2> legendre(unsigned n, geom::real x);

//------------IMPLEMENTATIONS---------------------------------------------------

template<class T>
constexpr T pow(const T& x, unsigned n)
{
  return (n == 0) ? T(1.) : x * pow(x, n - 1);
}

}

#endif // _LEGENDRE_HPP_
