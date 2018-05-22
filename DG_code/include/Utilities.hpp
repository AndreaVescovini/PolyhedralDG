/*!
    @file   Utilities.hpp
    @author Andrea Vescovini
    @brief  Here are defined some useful functions
*/

#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_

//! Some utilities
namespace Utilities
{

/*!
    @brief Integer power of a real number.

    This function implements recursively through templates the integer power of
    a real number.

    @param x        A real number.
    @param exponent The exponent.
*/
template<typename T>
constexpr T pow(const T& x, unsigned exponent)
{
  return (exponent == 0) ? T(1.0) : x * pow(x, exponent - 1);
}

} // namespace Utilities

#endif // _UTILITIES_HPP_
