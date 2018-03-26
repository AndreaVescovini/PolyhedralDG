#include "Operators.hpp"

namespace dgfem
{

geom::real Stiff::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i).dot(fe.getPhiDer(t, q, j));
}

geom::real Mass::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i) * fe.getPhi(t, q, j);
}

Eigen::Vector3d GradPhi::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i);
}

geom::real Phi::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i);
}

geom::real PhiI::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i);
}

geom::real PhiJ::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, j);
}

Eigen::Vector3d JumpInt::operator()(const FeFaceInt& fe, unsigned i, unsigned side, unsigned q) const
{
  return fe.getPhi(side, q, i) * (1 - 2 * side) * fe.getNormal() ;
}

}
