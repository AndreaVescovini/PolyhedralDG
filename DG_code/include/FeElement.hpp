#ifndef _FE_ELEMENT_HPP_
#define _FE_ELEMENT_HPP_

#include "PolyDG.hpp"
#include "Polyhedron.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <iostream>
#include <vector>

namespace PolyDG
{

class FeElement
{
public:
  using Element = Polyhedron;

  // Constructor that takes a geometrical polyhedron elem, the degree of polynomials,
  // the consequent number of degrees of freedom dof, the possible composition
  // of polynomilas of degree less or equal to degree into monomials basisComposition,
  // a quadrature rule tetraRule.
  FeElement(const Element& elem, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule3D& tetraRule);

  // Default copy-constructor.
  FeElement(const FeElement&) = default;

  // Default move-constructor.
  FeElement(FeElement&&) = default;

  // Function that returns the a const reference to the geomentrical element.
  inline const Element& getElem() const;

  // Function that returns the number of tetrhedra that compose elem_.
  inline SizeType getTetrahedraNo() const;

  // Function that return the absolute value of the determinant of the jacobian of the i-th tetrahedron.
  inline Real getAbsDetJac(SizeType i) const;

  // Function that returns the number of degrees of freedom that in 3D is
  // dof = (degree+1)*(degree+2)*(degree+3)/(3!)
  inline unsigned getDof() const;

  // Function that returns the number of quadrature points of tetraRule_.
  inline SizeType getQuadPointsNo() const;

  // Function that returns the p-th quadrature point in the tetrahedron t.
  inline Eigen::Vector3d getQuadPoint(SizeType t, SizeType p) const;

  // Function that returns the i-th quadrature weight.
  inline Real getWeight(SizeType i) const;

  // Function that returns the value of the basis function f, computed at the quadrature node p in
  // the tetrahedron t.
  inline Real getPhi(SizeType t, SizeType p, SizeType f) const;

  // Function that returns the vector of the gradient of the basis function f, computed at the
  // quadrature node p in the tetrahedron t.
  inline const Eigen::Vector3d& getPhiDer(SizeType t, SizeType p, SizeType f) const;

  // Prints all the computed values of the basis functions.
  void printBasis(std::ostream& out = std::cout) const;

  // Prints all the computed values of the gradient of the basis functions.
  void printBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeElement() = default;

private:

  // Reference to the element over which the finite element is computed
  const Element& elem_;

  // Number of degrees of freedom that in 3D is dof = (degree+1)*(degree+2)*(degree+3)/(3!)
  unsigned dof_;

  // Vector that contains all the possible degrees of the monomials that give
  // polynomials of degree less or equal to degree_
  const std::vector<std::array<unsigned, 3>>& basisComposition_;

  // Reference to the tetrahedral rule used for the computation of the basis functions
  const QuadRule3D& tetraRule_;

  // Vector containing all the values of the basis functions over this element
  std::vector<Real> phi_;

  // Vector containing all the values of the gradient of the basis functions over
  // this element
  std::vector<Eigen::Vector3d> phiDer_;

  // Auxiliary function that performs the computation of the basis functinos and
  // fills phi_ and phiDer_.
  void compute_basis();

  // Auxiliary function that, given the number of the tetrahedron t, quadrature
  // point p and basis function f, return the index in which the corresponding value
  // is stored in phi_ and phiDer_
  inline SizeType sub2ind(SizeType t, SizeType p, SizeType f) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const FeElement::Element& FeElement::getElem() const
{
  return elem_;
}

inline SizeType FeElement::getTetrahedraNo() const
{
  return elem_.getTetrahedraNo();
}

inline Real FeElement::getAbsDetJac(SizeType i) const
{
  return elem_.getTetra(i).getAbsDetJacobian();
}

inline unsigned FeElement::getDof() const
{
  return dof_;
}

inline SizeType FeElement::getQuadPointsNo() const
{
  return tetraRule_.getPointsNo();
}

inline Eigen::Vector3d FeElement::getQuadPoint(SizeType t, SizeType p) const
{
  return elem_.getTetra(t).getMap() * tetraRule_.getPoint(p);
}

inline Real FeElement::getWeight(SizeType i) const
{
  return tetraRule_.getWeight(i);
}

inline Real FeElement::getPhi(SizeType t, SizeType p, SizeType f) const
{
  return phi_[sub2ind(t, p, f)];
}

inline const Eigen::Vector3d& FeElement::getPhiDer(SizeType t, SizeType p, SizeType f) const
{
  return phiDer_[sub2ind(t, p, f)];
}

inline SizeType FeElement::sub2ind(SizeType t, SizeType p, SizeType f) const
{
  return f + dof_ * (p + t * tetraRule_.getPointsNo());
}

} // namespace PolyDG

#endif // _FE_ELEMENT_HPP_
