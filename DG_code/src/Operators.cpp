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

Eigen::Vector3d JumpPhi::operator()(const FeFaceInt& fe, unsigned i, int side, unsigned q) const
{
  return fe.getPhi(side, q, i) * (1.0 - 2.0 * side) * fe.getNormal();
}

Eigen::Vector3d JumpPhi::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return fe.getPhi(q, i) * fe.getNormal();
}

Eigen::Vector3d AverGradPhi::operator()(const FeFaceInt& fe, unsigned i, int side, unsigned q) const
{
  return 0.5 * fe.getPhiDer(side, q, i);
}

Eigen::Vector3d AverGradPhi::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return fe.getPhiDer(q, i);
}

PenaltyScaling::PenaltyScaling(geom::real sigma)
  : sigma_{sigma} {}

geom::real PenaltyScaling::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return sigma_ * fe.getPenaltyParam();
}

geom::real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return sigma_ * fe.getPenaltyParam();
}

} // namespace dgfem
