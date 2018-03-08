#include "FeFaceExt.hpp"

namespace dgfem
{

FeFaceExt::FeFaceExt(const theFaceExt& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule)
  : face_{face}, order_{order}, dofNo_{dofNo} {}

unsigned FeFaceExt::getOrder() const
{
  return order_;
}

unsigned FeFaceExt::getDofNo() const
{
  return dofNo_;
}

}
