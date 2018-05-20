/*!
    @file   FeFaceAbs.hpp
    @author Andrea Vescovini
    @brief  Abstract base class that defines the restriction over faces of finite elements
*/

#ifndef _FE_FACE_ABS_HPP_
#define _FE_FACE_ABS_HPP_

#include "FaceAbs.hpp"
#include "PolyDG.hpp"
#include "QuadRule.hpp"
#include "QuadRuleManager.hpp"

#include <Eigen/Core>

#include <array>
#include <iostream>
#include <vector>

namespace PolyDG
{

/*!
    @brief Abstract base class that defines the restriction over faces of finite elements

    This class is the base class that defines the restriction over faces of
    finite elements.
*/

class FeFaceAbs
{
public:
  /*!
      @brief Constructor

      @param face             A geometrical FaceAbs.
      @param dof              The degree of polynomials for the basis functions.
      @param basisComposition The composition of polynomials of degree less or
                              equal to degree into monomials.
      @param triaRule         A quadrature rule over triangles.
  */
  FeFaceAbs(const FaceAbs& face, unsigned dof, const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule2D& triaRule);

  //! Copy constructor
  FeFaceAbs(const FeFaceAbs&) = default;

  //! Move-constructor
  FeFaceAbs(FeFaceAbs&&) = default;

  /*!
      @brief Get the number of degrees of freedom

      This function returns the number of degrees of freedom in the FeElement,
      being in 3D it is (degree+1)*(degree+2)*(degree+3)/(3!).
  */
  inline unsigned getDof() const;

  //! Get the id number of the element to which the face belongs from the side @a "In"
  inline unsigned getElemIn() const;

  //! Get the number of quadrature points
  inline SizeType getQuadPointsNo() const;

  /*!
      @brief Get a quadrature weight

      This functions returns the p-th quadrature weight

      @param p The index of the quadrature weight required, it can be 0,..,getQuadPointsNo() - 1.
  */
  inline Real getWeight(SizeType p) const;

  /*!
      @brief Get a quadrature point

      This functions returns the p-th quadrature point.

      @param p The index of the quadrature point required, it can be 0,..,getQuadPointsNo() - 1.
  */
   inline Eigen::Vector3d getQuadPoint(SizeType p) const;

  /*!
      @brief Get the value of the penalty parameter for the face

      The penalty parameter is defined as:
      \f$ \max\limits_{\kappa \in \{\kappa^+, \kappa^-\}} \big\{ \frac{r^2}{h_\kappa}\big\}
      \f$ for internal faces and \f$ \frac{r^2}{h_\kappa} \f$ for external
      faces, where \f$ r \f$ is the degree of the FeSpace, \f$ h_{\kappa} \f$ is the
      diameter of the element \f$ \kappa \f$ and \f$ \kappa^+, \kappa^- \f$ are
      the two elements @a "In" and @a "Out" sharing the internal face.
  */
  inline Real getPenaltyParam() const;

  //! Get the the measure of the area doubled
  inline Real getAreaDoubled() const;

  /*!
      @brief  Get the normal vector
      @return @c Eigen::Vector3d containing the unitary normal vector, outward with
                 respect to the Tetrhedron @a "In".
  */
  inline const Eigen::Vector3d& getNormal() const;

  //! Print the values of the basis functions over the face
  virtual void printBasis(std::ostream& out) const = 0;

  //! Print the values of the gradient of the basis functions over the face
  virtual void printBasisDer(std::ostream& out) const = 0;

  //! Destructor
  virtual ~FeFaceAbs() = default;

protected:
  //! The geometrical face
  const FaceAbs& face_;

  //! Number of degrees of freedom that in 3D is dof = (degree+1)*(degree+2)*(degree+3)/(3!)
  unsigned dof_;

  //! Possible degrees of the monomials that multiplied togheter give polynomials of degree less or equal to degree_
  const std::vector<std::array<unsigned, 3>>& basisComposition_;

  //! Triangular quadrature rule used for the computation of the basis functions
  const QuadRule2D& triaRule_;

  //! Values of the all basis functions over the quadrature points of this face
  std::vector<Real> phi_;

  //! Values of the all gradients of the basis functions over the quadrature points of this face
  std::vector<Eigen::Vector3d> phiDer_;

  //! Penalty parameter
  Real penaltyParam_;

  //! Evaluate the basis functions and their gradient and fill phi_ and phDer_
  virtual void compute_basis() = 0;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeFaceAbs::getDof() const
{
  return dof_;
}

inline unsigned FeFaceAbs::getElemIn() const
{
  return face_.getTetIn().getPoly().getId();
}

inline SizeType FeFaceAbs::getQuadPointsNo() const
{
  return triaRule_.getPointsNo();
}

inline Real FeFaceAbs::getWeight(SizeType p) const
{
  return triaRule_.getWeight(p);
}

inline Eigen::Vector3d FeFaceAbs::getQuadPoint(SizeType p) const
{
  return face_.getTetIn().getMap() * (QuadRuleManager::instance().getFaceMap(face_.getFaceNoTetIn()) *
                                     triaRule_.getPoint(p).homogeneous());
}

inline Real FeFaceAbs::getPenaltyParam() const
{
  return penaltyParam_;
}

inline Real FeFaceAbs::getAreaDoubled() const
{
  return face_.getAreaDoubled();
}

inline const Eigen::Vector3d& FeFaceAbs::getNormal() const
{
  return face_.getNormal();
}

} // namespace PolyDG

#endif // _FE_FACE_ABS_HPP_
