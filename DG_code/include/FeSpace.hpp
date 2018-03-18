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

// Constructor that takes a mesh Th, the order of polynomials order and the required
// degrees of exactness for quadrature formulas
  FeSpace(TheMesh& Th, unsigned order, unsigned quad3DDegree, unsigned quad2DDegree);

// Constructor that takes a mesh Th and the order of polynomials order, the quadrature
// formulas are chosen to fit with order.
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
// Reference to the mesh over which the FeSpace is built
  const TheMesh& Th_;

// Order of polynomials
  unsigned order_;

// Number of degrees of freedom that in 3D is dofNo = (order+1)*(order+2)*(order+3)/(3!)
  unsigned dofNo_;

  std::vector<std::array<unsigned, 3>> basisComposition_;
  std::vector<FeElement> feElements_;
  std::vector<FeFaceInt> feFacesInt_;
  std::vector<FeFaceExt> feFacesExt_;
  const QuadRuleManager::Rule3D& tetraRule_;
  const QuadRuleManager::Rule2D& triaRule_;

// Auxiliary function that computes basisComposition_
  void integerComposition();

// Auxiliary function that fills feElements_, feFacesInt_ and feFacesExt_.
  void initialize();

};

}

#endif // _FE_SPACE_HPP_
