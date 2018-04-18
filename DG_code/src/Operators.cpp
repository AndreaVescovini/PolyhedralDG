#include "Operators.hpp"

namespace PolyDG
{

PenaltyScaling::PenaltyScaling(Real sigma)
  : sigma_{sigma} {}

Function::Function(const fun3real& fun)
  : fun_{fun} {}

} // namespace PolyDG
