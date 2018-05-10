#include "FeSpace.hpp"
#include "QuadRuleManager.hpp"

namespace PolyDG
{

FeSpace::FeSpace(Mesh& Th, unsigned order, unsigned quad3DDegree, unsigned quad2DDegree)
  : Th_{Th}, order_{order}, dof_{(order + 1) * (order + 2) * (order + 3) / 6},
    tetraRule_{QuadRuleManager::instance().getTetraRule(quad3DDegree)},
    triaRule_ {QuadRuleManager::instance().getTriaRule(quad2DDegree)}
  {
    integerComposition();
    initialize();
  }

// ho solo il termine di stiffness con le derivate e non il termine di massa,
// quindi posso abbassare l'ordine della quadratura nei tetraedri
FeSpace::FeSpace(Mesh& Th, unsigned order)
  : FeSpace(Th, order, 2*(order-1), 2*order) {}

void FeSpace::integerComposition()
{
  basisComposition_.reserve(dof_);

  int nx = order_;
  while(nx >= 0)
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
  feElements_.reserve(Th_.getPolyhedraNo());
  for(SizeType i = 0; i < Th_.getPolyhedraNo(); i++)
    feElements_.emplace_back(Th_.getPolyhedron(i), dof_, basisComposition_, tetraRule_);

  feFacesExt_.reserve(Th_.getFacesExtNo());
  for(SizeType i = 0; i < Th_.getFacesExtNo(); i++)
    feFacesExt_.emplace_back(Th_.getFaceExt(i), order_, dof_, basisComposition_, triaRule_);

  feFacesInt_.reserve(Th_.getFacesIntNo());
  for(SizeType i = 0; i < Th_.getFacesIntNo(); i++)
    feFacesInt_.emplace_back(Th_.getFaceInt(i), order_, dof_, basisComposition_, triaRule_);
}

void FeSpace::printElemBasis(std::ostream& out) const
{
  out << "-------- BASIS OVER ELEMENTS --------" << '\n';

  for(const FeElement& el : feElements_)
    el.printBasis(out);

  out << "-------------------------------------" << std::endl;
}

void FeSpace::printElemBasisDer(std::ostream& out) const
{
  out << "-------- BASIS DERIVATIVE OVER ELEMENTS --------" << '\n';

  for(const FeElement& el : feElements_)
    el.printBasisDer(out);

  out << "------------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasis(std::ostream& out) const
{
  out << "-------- BASIS OVER INTERNAL FACES --------" << '\n';

  for(const FeFaceInt& face : feFacesInt_)
    face.printBasis(out);

  out << "-------- BASIS OVER EXTERNAL FACES --------" << '\n';

  for(const FeFaceExt& face : feFacesExt_)
    face.printBasis(out);

  out << "-------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasisDer(std::ostream& out) const
{
  out << "-------- BASIS DERIVATIVE OVER INTERNAL FACES --------" << '\n';

  for(const FeFaceInt& face : feFacesInt_)
    face.printBasisDer(out);

  out << "-------- BASIS DERIVATIVE OVER EXTERNAL FACES --------" << '\n';

  for(const FeFaceExt& face : feFacesExt_)
    face.printBasisDer(out);

  out << "------------------------------------------------------" << std::endl;
}

} // namespace PolyDG
