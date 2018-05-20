/*!
    @file   FeFaceInt.hpp
    @author Andrea Vescovini
    @brief  Class that defines the restriction over internal faces of finite elements
*/

#ifndef _FE_FACE_INT_HPP_
#define _FE_FACE_INT_HPP_

#include "FaceInt.hpp"
#include "FeFaceAbs.hpp"
#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace PolyDG
{

/*!
    @brief Class that defines the restriction over external faces of finite elements

    This class inherits from FeFaceAbs and defines the restriction over internal
    faces of finite elements. The basis function and their gradient are evaluated
    at the 3D quadrature nodes that are provieded through a QuadRule2D and mapped
    to the face through QuadRuleManager::getFaceMap
*/

class FeFaceInt : public FeFaceAbs
{
public:
  /*!
      @brief Constructor

      This constructor creates the FeFaceInt, computes the penalty parameter and
      computes the values of the basis functions and their gradient at the
      quadrature points.

      @param face             A geometrical FaceInt.
      @param degree           The degree of the FeSpace.
      @param dof              The degree of polynomials for the basis functions.
      @param basisComposition The composition of polynomials of degree less or
                              equal to degree into monomials.
      @param triaRule         A quadrature rule over triangles.
  */
  FeFaceInt(const FaceInt& face, unsigned degree, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule2D& triaRule);

  //! Copy constructor
  FeFaceInt(const FeFaceInt&) = default;

  //! Move constructor
  FeFaceInt(FeFaceInt&&) = default;

  /*!
      @brief Get the value of the basis function

      This functions returns the value of the f-th basis function at the
      p-th quadrature, from the side s.

      @param s The SideType of from which the function is taken.
      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline Real getPhi(SideType s, SizeType p, SizeType f) const;

  /*!
      @brief Get the value of the gradient of the basis function

      This functions returns the value of the gradient of the f-th basis function
      at the p-th quadrature, from the side s

      @param s The SideType of from which the function is taken.
      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline const Eigen::Vector3d& getPhiDer(SideType s, SizeType p, SizeType f) const;

  //! Get the id number of the element to which the face belongs from the side @a "Out"
  inline unsigned getElemOut() const;

  //! Print all the computed values of the basis functions, from both sides
  void printBasis(std::ostream& out) const override;

  //! Print all the computed values of the gradient of the basis functions from both sides
  void printBasisDer(std::ostream& out) const override;

  //! Destructor
  virtual ~FeFaceInt() = default;

private:
  //! Evaluate the basis functions and their gradient at the quadrature nodes and fill phi_ and phiDer_
  void compute_basis() override;

  //! Combine the three indices of a side s, quadrature point and a basis function into one index
  inline SizeType sub2ind(SideType s, SizeType p, SizeType f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceInt::getPhi(SideType s, SizeType p, SizeType f) const
{
  return phi_[sub2ind(s, p, f)];
}

inline const Eigen::Vector3d& FeFaceInt::getPhiDer(SideType s, SizeType p, SizeType f) const
{
  return phiDer_[sub2ind(s, p, f)];
}

inline unsigned FeFaceInt::getElemOut() const
{
  return static_cast<const FaceInt&>(face_).getTetOut().getPoly().getId();
}

inline SizeType FeFaceInt::sub2ind(SideType s, SizeType p, SizeType f) const
{
  return (s == In) + 2 * (f + p * dof_);
}

} // namespace PolyDG

#endif // _FE_FACE_INT_HPP_
