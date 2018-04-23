#ifndef _LEGENDRE_BASIS_HPP_
#define _LEGENDRE_BASIS_HPP_

#include "PolyDG.hpp"

#include <iostream>

namespace PolyDG
{

// Function that return the value of the legendre polynomial of order n at the
// point x belonging to [-1, 1] and the value of the derivative at the same point.

template<unsigned N>
Real LegendreBasis(Real x);

template<unsigned N>
Real LegendreBasisDer(Real x);

// Recursive implementation through templates of the integer power of a real number
template<typename T>
constexpr T pow(const T& x, unsigned n);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template<unsigned N>
Real LegendreBasis(Real /*x*/)
{
  std::cerr << "The basis function of order " << N << " are not implemented" << std::endl;
  return 0.0;
}

template<unsigned N>
Real LegendreBasisDer(Real /*x*/)
{
  std::cerr << "The basis function of order " << N << " are not implemented" << std::endl;
  return 0.0;
}

template<typename T>
constexpr T pow(const T& x, unsigned n)
{
  return (n == 0) ? T(1.0) : x * pow(x, n - 1);
}

} // namespace PolyDG

#endif // _LEGENDRE_BASIS_HPP_
