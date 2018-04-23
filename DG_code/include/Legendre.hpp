#ifndef _LEGENDRE_HPP_
#define _LEGENDRE_HPP_

#include "PolyDG.hpp"

#include <array>

namespace basis
{

// Recursive implementation through templates of the integer power of a real number
template<typename T>
constexpr T pow(const T& x, unsigned n);

// Function that return the value of the legendre polynomial of order n at the
// point x belonging to [-1, 1] and the value of the derivative at the same point.
std::array<PolyDG::Real, 2> legendre(unsigned n, PolyDG::Real x);


//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template<typename T>
constexpr T pow(const T& x, unsigned n)
{
  return (n == 0) ? T(1.) : x * pow(x, n - 1);
}

}

#endif // _LEGENDRE_HPP_
