/*!
    @file   FeFaceExt.hpp
    @author Andrea Vescovini
    @brief  Class that defines the restriction over external faces of finite elements
*/

#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FaceExt.hpp"
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

    This class inherits from FeFaceAbs and defines the restriction over external
    faces of finite elements. The basis function and their gradient are evaluated
    at the 3D quadrature nodes that are provieded through a QuadRule2D and mapped
    to the face through QuadRuleManager::getFaceMap
*/

class FeFaceExt : public FeFaceAbs
{
public:
  /*!
      @brief Constructor

      This constructor creates the FeFaceExt, computes the penalty parameter and
      computes the values of the basis functions and their gradient at the
      quadrature points.

      @param face             A geometrical FaceExt.
      @param degree           The degree of the FeSpace.
      @param dof              The degree of polynomials for the basis functions.
      @param basisComposition The composition of polynomials of degree less or
                              equal to degree into monomials.
      @param triaRule         A quadrature rule over triangles.
  */
  FeFaceExt(const FaceExt& face, unsigned degree, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule2D& triaRule);

  //! Copy constructor
  FeFaceExt(const FeFaceExt&) = default;

  //! Move constructor
  FeFaceExt(FeFaceExt&&) = default;

  /*!
      @brief Get the value of the basis function

      This functions returns the value of the f-th basis function at the
      p-th quadrature.

      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline Real getPhi(SizeType p, SizeType f) const;

  /*!
      @brief Get the value of the gradient of the basis function

      This functions returns the value of the gradient of the f-th basis function
      at the p-th quadrature.

      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
      @param f The index of the basis function required, it can be 0,...,getDof() - 1.
  */
  inline const Eigen::Vector3d& getPhiDer(SizeType p, SizeType f) const;

  //! Get the label related to the FaceExt
  inline BCType getBClabel() const;

  //! Print all the computed values of the basis functions
  void printBasis(std::ostream& out) const override;

  //! Print all the computed values of the gradient of the basis functions
  void printBasisDer(std::ostream& out) const override;

  //! Destructor
  virtual ~FeFaceExt() = default;

private:
  //! Evaluate the basis functions and their gradient at the quadrature nodes and fill phi_ and phiDer_
  void compute_basis() override;

  //! Combine the two indices of a quadrature point and a basis function into one index
  inline SizeType sub2ind(SizeType p, SizeType f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceExt::getPhi(SizeType p, SizeType f) const
{
  return phi_[sub2ind(p, f)];
}

inline const Eigen::Vector3d& FeFaceExt::getPhiDer(SizeType p, SizeType f) const
{
  return phiDer_[sub2ind(p, f)];
}

inline BCType FeFaceExt::getBClabel() const
{
  return static_cast<const FaceExt&>(face_).getBClabel();
}

inline SizeType FeFaceExt::sub2ind(SizeType p, SizeType f) const
{
  return f + dof_ * p;
}

} // namespace PolyDG

#endif // _FE_FACE_EXT_HPP_
