#include "FeSpace.hpp"

namespace dgfem {

FeSpace::FeSpace(TheMesh& Th, unsigned order, unsigned quad3DDegree, unsigned quad2DDegree)
  : Th_{Th}, order_{order},
    tetraRule_{QuadRuleManager::getTetraRule(quad3DDegree)},
    triaRule_ {QuadRuleManager::getTriaRule(quad2DDegree)}
  {
    dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;
    integerComposition();

    initialize();
  }

// ho solo il termine di stiffness con le derivate e non il termine di massa,
// quindi posso abbassare l'ordine della quadratura nei tetraedri
FeSpace::FeSpace(TheMesh& Th, unsigned order)
  : FeSpace(Th, order, 2*(order-1), 2*order) {}

void FeSpace::setOrder(unsigned order)
{
  order_ = order;
  dofNo_ = (order + 1) * (order + 2) * (order + 3) / 6.;
  integerComposition();
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
    feElements_.emplace_back(Th_.getPolyhedron(i), dofNo_,
                             basisComposition_, tetraRule_);
  }

  elemNo = Th_.getFacesExtNo();
  feFacesExt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesExt_.emplace_back( Th_.getFaceExt(i), order_, dofNo_, basisComposition_, triaRule_);
  }

  elemNo = Th_.getFacesIntNo();
  feFacesInt_.reserve(elemNo);
  for(unsigned i = 0; i < elemNo; i++)
  {
    feFacesInt_.emplace_back(Th_.getFaceInt(i), order_, dofNo_, basisComposition_, triaRule_);
  }
}

void FeSpace::printElemBasis(std::ostream& out) const
{
  out << "-------- BASIS OVER ELEMENTS --------" << "\n";

  for(auto& el : feElements_)
  {
    el.printBasis(out);
  }

  out << "-------------------------------------" << std::endl;
}

void FeSpace::printElemBasisDer(std::ostream& out) const
{
  out << "-------- BASIS DERIVATIVE OVER ELEMENTS --------" << "\n";

  for(auto& el : feElements_)
  {
    el.printBasisDer(out);
  }

  out << "------------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasis(std::ostream& out) const
{
  out << "-------- BASIS OVER INTERNAL FACES --------" << "\n";

  for(auto& face : feFacesInt_)
  {
    face.printBasis(out);
  }

  out << "-------- BASIS OVER EXTERNAL FACES --------" << "\n";

  for(auto& face : feFacesExt_)
  {
    face.printBasis(out);
  }

  out << "-------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasisDer(std::ostream& out) const
{
  out << "-------- BASIS DERIVATIVE OVER INTERNAL FACES --------" << "\n";

  for(auto& face : feFacesInt_)
  {
    face.printBasisDer(out);
  }

  out << "-------- BASIS DERIVATIVE OVEREXTERNAL FACES --------" << "\n";

  for(auto& face : feFacesExt_)
  {
    face.printBasisDer(out);
  }

  out << "------------------------------------------------------" << std::endl;
}

}
