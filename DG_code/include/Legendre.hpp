/*!
    @file   Legendre.hpp
    @author Andrea Vescovini
    @brief  Here are defined Legendre polynomials
*/

#ifndef _LEGENDRE_HPP_
#define _LEGENDRE_HPP_

#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Evaulate a Legendre polynomial

    This function evaulates the Legendre polynomial of degree n at the point x
    belonging to [-1, 1].

    @warning It is implemented only for n <= 8, for n > 8 a @c std::domain_error
             exception is thrown.
    @param n Degree of the polynomial.
    @param x Evaluation point, it must be in [-1, 1].
*/
Real legendre(unsigned n, PolyDG::Real x);

/*!
    @brief Evaulate a the first derivative of a Legendre polynomial

    This function evaulates the first derivative of the Legendre polynomial of
    degree n at the point x belonging to [-1, 1].

    @warning It is implemented only for n <= 8, for n > 8 a @c std::domain_error
             exception is thrown.
    @param n Degree of the polynomial.
    @param x Evaluation point, it must be in [-1, 1].
*/
Real legendreDer(unsigned n, PolyDG::Real x);

/*!
    @brief Integer power of a real number.

    This function implements recursively through templates the integer power of
    a real number.

    @param x        A real number.
    @param exponent The exponent.
*/
template<typename T>
constexpr T pow(const T& x, unsigned exponent);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template<typename T>
constexpr T pow(const T& x, unsigned exponent)
{
  return (exponent == 0) ? T(1.0) : x * pow(x, exponent - 1);
}

}

#endif // _LEGENDRE_HPP_
