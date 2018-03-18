#ifndef _FE_ELEMENT_HPP_
#define _FE_ELEMENT_HPP_

#include "Polyhedron.hpp"
#include "geom.hpp"
#include "QuadRuleManager.hpp"
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <iostream>

namespace dgfem
{

class FeElement
{
public:
  using TheElem = geom::Polyhedron;

  FeElement(const TheElem& elem, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule3D& tetraRule);

// Returns the order of the polynomials emplyed in the finite element
  unsigned getOrder() const;

// Returns the number of degrees of freedom that in 3D is
// dofNo = (order+1)*(order+2)*(order+3)/(3!)
  unsigned getDofNo() const;

// Returns the value of the basis function f, computed at the quadrature node p in
// the tetrahedron t
  geom::real getPhi(unsigned t, unsigned p, unsigned f) const;

// Returns the vector of the gradient of the basis function f, computed at the
// quadrature node p in the tetrahedron t
  const Eigen::Vector3d& getPhiDer(unsigned t, unsigned p, unsigned f) const;

// Prints all the computed values of the basis functions
  void printBasis(std::ostream& out = std::cout) const;

// Prints all the computed values of the gradient of the basis functions
  void printBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeElement() = default;

private:

// Reference to the element over which the finite element is computed 
  const TheElem& elem_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
  const std::vector<std::array<unsigned, 3>>& basisComposition_;
  const QuadRuleManager::Rule3D& tetraRule_;

  std::vector<geom::real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

  void compute_basis();
  unsigned sub2ind(unsigned t, unsigned p, unsigned f) const;
};

}

#endif // _FE_ELEMENT_HPP_
