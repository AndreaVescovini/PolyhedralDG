#ifndef _FE_ELEMENT_HPP_
#define _FE_ELEMENT_HPP_

#include "Polyhedron.hpp"
#include "geom.hpp"
#include "QuadRuleManager.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <iostream>

namespace dgfem
{

class FeElement
{
public:
  using TheElem = geom::Polyhedron;

// Constructor that takes a geometrical polyhedron elem, the order of polynomials
// order, the consequent number of degrees of freedom dofNo, the possible composition
// of polynomilas of degree less or equal to order into monomials basisComposition,
// a quadrature rule tetraRule.
  FeElement(const TheElem& elem, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule3D& tetraRule);

// Returns the number of degrees of freedom that in 3D is
// dofNo = (order+1)*(order+2)*(order+3)/(3!)
  unsigned getDofNo() const;

// Returns the number of tetrhedra that compose elem_
  unsigned getTetrahedraNo() const;

// Returns the number of quadrature points of tetraRule_
  unsigned getQuadPointsNo() const;

// Return the absolute value of the determinant of the jacobian of the i-th tetrahedron
  geom::real getAbsDetJac(unsigned i) const;

// Returns the value of the basis function f, computed at the quadrature node p in
// the tetrahedron t
  geom::real getPhi(unsigned t, unsigned p, unsigned f) const;

// Returns the vector of the gradient of the basis function f, computed at the
// quadrature node p in the tetrahedron t
  const Eigen::Vector3d& getPhiDer(unsigned t, unsigned p, unsigned f) const;

// Returns the i-th quadrature weight
  geom::real getWeight(unsigned i) const;

// Prints all the computed values of the basis functions
  void printBasis(std::ostream& out = std::cout) const;

// Prints all the computed values of the gradient of the basis functions
  void printBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeElement() = default;

private:

// Reference to the element over which the finite element is computed
  const TheElem& elem_;

// Number of degrees of freedom that in 3D is dofNo = (order+1)*(order+2)*(order+3)/(3!)
  unsigned dofNo_;

// Vector that contains all the possible degrees of the monomials that give
// polynomials of degree less or equal to order_
  const std::vector<std::array<unsigned, 3>>& basisComposition_;

// Reference to the tetrahedral rule used for the computation of the basis functions
  const QuadRuleManager::Rule3D& tetraRule_;

// Vector containing all the values of the basis functions over this element
  std::vector<geom::real> phi_;

// Vector containing all the values of the gradient of the basis functions over
// this element
  std::vector<Eigen::Vector3d> phiDer_;

// Auxiliary function that performs the computation of the basis functinos and
// fills phi_ and phiDer_
  void compute_basis();

// Auxiliary function that, given the number of the tetrahedron t, quadrature
// point p and basis function f, return the index in which the corresponding value
// is stored in phi_ and phiDer_
  unsigned sub2ind(unsigned t, unsigned p, unsigned f) const;
};

}

#endif // _FE_ELEMENT_HPP_
