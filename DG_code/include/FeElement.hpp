/*!
    @file   FeElement.hpp
    @author Andrea Vescovini
    @brief  Class that defines finite elements over polyhedra
*/

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

/*!
    @brief Class that defines finite elements over polyhedra

    This class defines finite elements over polyhedra. The basis function and
    their gradient are evaluated at the 3D quadrature nodes that are provieded
    through a QuadRule3D in the constructor.
*/

class FeElement
{
public:
  //! Alias for Polyhedron
  using Element = Polyhedron;

  /*!
      @brief Constructor

      This constructor creates the FeElement and computes the values of the
      basis functions and their gradient at the quadrature points.

      @param elem             A geometrical Polyhedron.
      @param dof              The degree of polynomials for the basis functions.
      @param basisComposition The composition of polynomials of degree less or
                              equal to degree into monomials.
      @param tetraRule        A quadrature rule over tetrahedra.
  */
  FeElement(const Element& elem, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule3D& tetraRule);

  //! Copy constructor
  FeElement(const FeElement&) = default;

  //! Move constructor
  FeElement(FeElement&&) = default;

  //! Get the geometrical element
  inline const Element& getElem() const;

  //! Get the number of tetrhedra that compose the Element
  inline SizeType getTetrahedraNo() const;

  /*!
      @brief Get the absolute value of the determinant of a jacobian

      This functions returns the absolute value of the determinant of the
      jacobian of the i-th Tetrahedron.

      @param i The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1.
  */
  inline Real getAbsDetJac(SizeType i) const;

  /*!
      @brief Get the number of degrees of freedom

      This function returns the number of degrees of freedom in the FeElement,
      being in 3D it is (degree+1)*(degree+2)*(degree+3)/(3!).
  */
  inline unsigned getDof() const;

  //! Get the number of quadrature points
  inline SizeType getQuadPointsNo() const;

  /*!
      @brief Get a quadrature point

      This functions returns the p-th quadrature point in the tetrahedron t.

      @param t The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1.
      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d getQuadPoint(SizeType t, SizeType p) const;

  /*!
      @brief Get a quadrature weight

      This functions returns the i-th quadrature weight

      @param i The index of the quadrature weight required, it can be 0,..,getQuadPointsNo() - 1.
  */
  inline Real getWeight(SizeType i) const;

  /*!
      @brief Get the value of the basis function

      This functions returns the value of the f-th basis function at the
      p-th quadrature point in the tetrahedron t.

      @param t The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1.
      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline Real getPhi(SizeType t, SizeType p, SizeType f) const;

  /*!
      @brief Get the value of the gradient of the basis function

      This functions returns the value of the gradient of the f-th basis function
      at the p-th quadrature point in the tetrahedron t.

      @param t The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1.
      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline const Eigen::Vector3d& getPhiDer(SizeType t, SizeType p, SizeType f) const;

  //! Prints all the computed values of the basis functions
  void printBasis(std::ostream& out = std::cout) const;

  //! Prints all the computed values of the gradient of the basis functions
  void printBasisDer(std::ostream& out = std::cout) const;

  //! Destructor
  virtual ~FeElement() = default;

private:

  //! Element over which the finite element is computed
  const Element& elem_;

  //! Number of degrees of freedom that in 3D is dof = (degree+1)*(degree+2)*(degree+3)/(3!)
  unsigned dof_;

  //! Possible degrees of the monomials that multiplied togheter give polynomials of degree less or equal to degree_
  const std::vector<std::array<unsigned, 3>>& basisComposition_;

  //! Tetrahedral quadrature rule used for the computation of the basis functions
  const QuadRule3D& tetraRule_;

  //! Values of the all basis functions over the quadrature points of this Element
  std::vector<Real> phi_;

  //! Values of the all gradients of the basis functions over the quadrature points of this Element
  std::vector<Eigen::Vector3d> phiDer_;

  //! Evaluate the basis functions and their gradient at the quadrature nodes and fill phi_ and phiDer_
  void compute_basis();

  //! Combine the three indices of a Tetrahedron, a quadrature point and a basis function into one index
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
