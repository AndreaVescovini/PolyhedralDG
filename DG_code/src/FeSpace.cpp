#include "FeSpace.hpp"
#include "QuadRuleManager.hpp"

namespace dgfem {

 // ho solo il termine di stiffness con le derivate e non il termine di massa,
 // quindi posso abbassare l'ordine della quadratura nei tetraedri
FeSpace::FeSpace(theMesh& Th, unsigned order)
  : Th_{Th}, order_{order},
    tetraRule_{QuadRuleManager::getTetraRule(2* (order - 1))},
    triaRule_ {QuadRuleManager::getTriaRule(2 * order)}
  {
    dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;
    integerComposition();



    initialize();
  }

void FeSpace::setOrder(unsigned order)
{
  order_ = order;
  dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;
  integerComposition();
}

unsigned FeSpace::getOrder() const
{
  return order_;
}

unsigned FeSpace::getDofNo() const
{
  return dofNo_;
}

void FeSpace::integerComposition()
{
  basisComposition_.reserve(dofNo_);

  int nx = order_;
  while (nx >= 0)
  {
    int ny = order_ - nx;
    while(ny >= 0)
    {
      int nz = order_ - nx - ny;
      while(nz >= 0)
      {
        basisComposition_.emplace_back(std::array<unsigned, 3>{static_cast<unsigned>(nx),
                                                               static_cast<unsigned>(ny),
                                                               static_cast<unsigned>(nz)});
        nz--;
      }
      ny--;
    }
    nx--;
  }
}

void FeSpace::initialize()
{
  unsigned elemNo = Th_.getPolyhedraNo();
  feElements_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feElements_.emplace_back(Th_.getPolyhedron(i), order_, dofNo_,
                             basisComposition_, tetraRule_);
  }

  elemNo = Th_.getFacesExtNo();
  feFacesExt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesExt_.emplace_back(Th_.getFaceExt(i), order_, dofNo_,
                             basisComposition_, triaRule_);
  }

  elemNo = Th_.getFacesIntNo();
  feFacesInt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesInt_.emplace_back(Th_.getFaceInt(i), order_, dofNo_,
                             basisComposition_, triaRule_);
  }
}

}
