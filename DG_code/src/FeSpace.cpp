#include "FeSpace.hpp"

namespace dgfem {

FeSpace::FeSpace(theMesh& Th, unsigned order)
  : Th_{Th}, order_{order}
  {
    dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;

    initialize();
  }

void FeSpace::setOrder(unsigned order)
{
  order_ = order;
  dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;
}

unsigned FeSpace::getOrder() const
{
  return order_;
}

unsigned FeSpace::getDofNo() const
{
  return dofNo_;
}

void FeSpace::initialize()
{
  unsigned elemNo = Th_.getPolyhedraNo();
  feElements_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feElements_.emplace_back(Th_.getPolyhedron(i));
  }

  elemNo = Th_.getFacesExtNo();
  feFacesExt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesExt_.emplace_back(Th_.getFaceExt(i));
  }

  elemNo = Th_.getFacesIntNo();
  feFacesInt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesInt_.emplace_back(Th_.getFaceInt(i));
  }
}

}
