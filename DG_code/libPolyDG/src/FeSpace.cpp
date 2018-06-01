/*!
    @file   FeSpace.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class FeSpace
*/

#include "FeSpace.hpp"
#include "QuadRuleManager.hpp"
#include "Watch.hpp"

namespace PolyDG
{

FeSpace::FeSpace(Mesh& Th, unsigned degree, unsigned doeQuad3D, unsigned doeQuad2D)
  : Th_{Th}, degree_{degree}, dof_{(degree + 1) * (degree + 2) * (degree + 3) / 6},
    tetraRule_{QuadRuleManager::instance().getTetraRule(doeQuad3D)},
    triaRule_ {QuadRuleManager::instance().getTriaRule(doeQuad2D)}
  {
    integerComposition();
    initialize();
  }

// ho solo il termine di stiffness con le derivate e non il termine di massa,
// quindi posso abbassare l'ordine della quadratura nei tetraedri
FeSpace::FeSpace(Mesh& Th, unsigned degree)
  : FeSpace(Th, degree, 2 * (degree - 1), 2 * degree) {}

void FeSpace::integerComposition()
{
  basisComposition_.reserve(dof_);

  int nx = degree_;
  while(nx >= 0)
  {
    int ny = degree_ - nx;
    while(ny >= 0)
    {
      int nz = degree_ - nx - ny;
      while(nz >= 0)
      {
        basisComposition_.emplace_back(std::array<unsigned, 3>{{static_cast<unsigned>(nx),
                                                                static_cast<unsigned>(ny),
                                                                static_cast<unsigned>(nz)}});
        nz--;
      }
      ny--;
    }
    nx--;
  }
}

void FeSpace::initialize()
{
  #ifdef VERBOSITY
    std::cout << "Inizializing FeElements...";
    Utilities::Watch ch;
    ch.start();
  #endif

  feElements_.reserve(Th_.getPolyhedraNo());
  for(SizeType i = 0; i < Th_.getPolyhedraNo(); i++)
    feElements_.emplace_back(Th_.getPolyhedron(i), dof_, basisComposition_, tetraRule_);

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << "\nInizializing FeFacesExt...";
    ch.reset();
    ch.start();
  #endif

  feFacesExt_.reserve(Th_.getFacesExtNo());
  for(SizeType i = 0; i < Th_.getFacesExtNo(); i++)
    feFacesExt_.emplace_back(Th_.getFaceExt(i), degree_, dof_, basisComposition_, triaRule_);

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << "\nInizializing FeFacesInt...";
    ch.reset();
    ch.start();
  #endif

  feFacesInt_.reserve(Th_.getFacesIntNo());
  for(SizeType i = 0; i < Th_.getFacesIntNo(); i++)
    feFacesInt_.emplace_back(Th_.getFaceInt(i), degree_, dof_, basisComposition_, triaRule_);

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
}

void FeSpace::printInfo(std::ostream& out) const
{
  out << "-------------------- FESPACE INFO --------------------" << '\n';
  out << "Degree = " << degree_ << '\n';
  out << "Degrees of freedom per element: " << dof_ << '\n';
  out << "Elements: " << feElements_.size() << '\n';
  out << "Total degrees of freedom: " << dof_ * feElements_.size() <<'\n';
  out << "Quadrature Rule 3D: degree of exactness = " << tetraRule_.getDoe() << ", points: " << tetraRule_.getPointsNo() << '\n';
  out << "Quadrature Rule 2D: degree of exactness = " << triaRule_.getDoe() << ", points: " << triaRule_.getPointsNo() << '\n';
  out << "------------------------------------------------------" << std::endl;
}

void FeSpace::printElemBasis(std::ostream& out) const
{
  out << "---------------- BASIS OVER ELEMENTS -----------------" << '\n';

  for(const FeElement& el : feElements_)
    el.printBasis(out);

  out << "------------------------------------------------------" << std::endl;
}

void FeSpace::printElemBasisDer(std::ostream& out) const
{
  out << "----------- BASIS DERIVATIVE OVER ELEMENTS -----------" << '\n';

  for(const FeElement& el : feElements_)
    el.printBasisDer(out);

  out << "------------------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasis(std::ostream& out) const
{
  out << "------------- BASIS OVER EXTERNAL FACES --------------" << '\n';

  for(const FeFaceExt& face : feFacesExt_)
    face.printBasis(out);

  out << "------------- BASIS OVER INTERNAL FACES --------------" << '\n';

  for(const FeFaceInt& face : feFacesInt_)
    face.printBasis(out);

  out << "------------------------------------------------------" << std::endl;
}

void FeSpace::printFaceBasisDer(std::ostream& out) const
{
  out << "-------- BASIS DERIVATIVE OVER EXTERNAL FACES --------" << '\n';

  for(const FeFaceExt& face : feFacesExt_)
    face.printBasisDer(out);

  out << "-------- BASIS DERIVATIVE OVER INTERNAL FACES --------" << '\n';

  for(const FeFaceInt& face : feFacesInt_)
    face.printBasisDer(out);

  out << "------------------------------------------------------" << std::endl;
}

} // namespace PolyDG
