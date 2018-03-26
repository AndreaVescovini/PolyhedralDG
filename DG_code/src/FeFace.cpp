#include "FeFace.hpp"

namespace dgfem
{

FeFace::FeFace(unsigned order, unsigned dofNo,
               const std::vector<std::array<unsigned, 3>>& basisComposition,
               const QuadRule<Eigen::Vector2d>& triaRule)
  : dofNo_{dofNo}, basisComposition_{basisComposition}, triaRule_{triaRule} {}

unsigned FeFace::getDofNo() const
{
  return dofNo_;
}

unsigned FeFace::getQuadPointsNo() const
{
  return triaRule_.getPointsNo();
}

geom::real FeFace::getWeight(unsigned i) const
{
  return triaRule_.getWeight(i);
}

geom::real FeFace::getPenaltyParam() const
{
  return penaltyParam_;
}

}
