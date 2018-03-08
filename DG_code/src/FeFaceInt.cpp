#include "FeFaceInt.hpp"

namespace dgfem
{

FeFaceInt::FeFaceInt(const theFaceInt& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule)
  : face_{face}, order_{order}, dofNo_{dofNo} {}

unsigned FeFaceInt::getOrder() const
{
  return order_;
}

unsigned FeFaceInt::getDofNo() const
{
  return dofNo_;
}

}
