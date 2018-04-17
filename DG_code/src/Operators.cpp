#include "Operators.hpp"

namespace dgfem
{

PenaltyScaling::PenaltyScaling(geom::real sigma)
  : sigma_{sigma} {}

Function::Function(const fun3real& fun)
  : fun_{fun} {}

} // namespace dgfem
