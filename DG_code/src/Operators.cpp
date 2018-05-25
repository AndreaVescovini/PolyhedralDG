/*!
    @file   Operators.cpp
    @author Andrea Vescovini
    @brief  Implementation of the classes for the operators
*/

#include "Operators.hpp"

namespace PolyDG
{

PenaltyScaling::PenaltyScaling(Real sigma)
  : sigma_{sigma} {}

Function::Function(const funR3R1& fun)
  : fun_{fun} {}

Function3::Function3(const funR3R3& fun)
  : fun_{fun} {}

} // namespace PolyDG
