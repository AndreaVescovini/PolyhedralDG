#ifndef _FE_SPACE_HPP_
#define _FE_SPACE_HPP_

#include "Mesh.hpp"
#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "QuadRuleManager.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <iostream>

namespace dgfem {

class FeSpace
{
public:
  using TheMesh = geom::Mesh;

  FeSpace(TheMesh& Th, unsigned order, unsigned quad3DDegree, unsigned quad2DDegree);
  FeSpace(TheMesh& Th, unsigned order);

  void setOrder(unsigned order);
  unsigned getOrder() const;
  unsigned getDofNo() const;

  void printElemBasis(std::ostream& out = std::cout) const;
  void printElemBasisDer(std::ostream& out = std::cout) const;
  void printFaceBasis(std::ostream& out = std::cout) const;
  void printFaceBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeSpace() = default;

private:
  const TheMesh& Th_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
  std::vector<std::array<unsigned, 3>> basisComposition_;
  std::vector<FeElement> feElements_;
  std::vector<FeFaceInt> feFacesInt_;
  std::vector<FeFaceExt> feFacesExt_;
  const QuadRuleManager::Rule3D& tetraRule_;
  const QuadRuleManager::Rule2D& triaRule_;

  void integerComposition();
  void initialize();


};

}

#endif // _FE_SPACE_HPP_
