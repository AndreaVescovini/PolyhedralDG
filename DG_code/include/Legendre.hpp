/*!
    @file   Legendre.hpp
    @author Andrea Vescovini
    @brief  Here are defined Legendre polynomials and their derivative
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

} // namespace PolyDG

#endif // _LEGENDRE_HPP_
