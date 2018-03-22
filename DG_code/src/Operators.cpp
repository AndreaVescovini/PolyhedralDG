#include "Operators.hpp"

namespace dgfem
{

geom::real Stiff::operator()(const FeElement& fe, unsigned i, unsigned j,
                             unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, j).transpose() * fe.getPhiDer(t, q, i);
}

}
