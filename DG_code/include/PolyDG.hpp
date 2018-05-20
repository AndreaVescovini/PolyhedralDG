/*!
    @file   PolyDG.hpp
    @author Andrea Vescovini
    @brief  Here some alias and enums are defined
*/

#ifndef _POLY_DG_HPP_
#define _POLY_DG_HPP_

#include <cstddef>

//! Library namespace
namespace PolyDG
{

//! Alias for real numbers
using Real = double;

//! Alias for size type of stl containers
using SizeType = std::size_t;

//! Alias for the type of boundary condition
using BCType = int;

//! Enum for the two sides of faces.
enum SideType { Out, In };

} // namespace PolyDG

#endif // _POLY_DG_HPP_
