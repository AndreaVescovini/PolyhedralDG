#ifndef _LEGENDRE_HPP_
#define _LEGENDRE_HPP_

#include "PolyDG.hpp"

namespace PolyDG
{

// Function that return the value of the Legendre polynomial of order n at the
// point x belonging to [-1, 1]. It is implemented only for n<=6.
Real legendre(unsigned n, PolyDG::Real x);

// Function that return the value of the derivative of the Legendre polynomial
// of order n at the point x belonging to [-1, 1]. It is implemented only for n<=6.
Real legendreDer(unsigned n, PolyDG::Real x);

// Recursive implementation through templates of the integer power of a real number.
template<typename T>
constexpr T pow(const T& x, unsigned n);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template<typename T>
constexpr T pow(const T& x, unsigned n)
{
  return (n == 0) ? T(1.0) : x * pow(x, n - 1);
}

}

#endif // _LEGENDRE_HPP_
