#ifndef _FE_FACE_HPP_
#define _FE_FACE_HPP_

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

class FeFace
{
public:
  FeFace(const FaceAbs& face, unsigned dof, const std::vector<std::array<unsigned, 3>>& basisComposition,
         const QuadRule2D& triaRule);

  // Default copy-constructor.
  FeFace(const FeFace&) = default;

  // Default move-constructor.
  FeFace(FeFace&&) = default;

  // Function that returns the number of degrees of freedom.
  inline unsigned getDof() const;

  // Function that returns the id-number of the elements to which the face belongs
  // from inside.
  inline unsigned getElemIn() const;

  // Function that returns the number of quadrature points currently used over
  // the face.
  inline SizeType getQuadPointsNo() const;

  // Function that returns the p-th quadrature weight.
  inline Real getWeight(SizeType p) const;

  // Function returning the p-th quarature point.
  inline Eigen::Vector3d getQuadPoint(SizeType p) const;

  // Function returning the value of the penalty parameter for the face.
  inline Real getPenaltyParam() const;

  // Function returning the doubled area of the face.
  inline Real getAreaDoubled() const;

  // Functino returning the normal to the face.
  inline const Eigen::Vector3d& getNormal() const;

  // Function that prints the values of the basis functions over the face.
  virtual void printBasis(std::ostream& out) const = 0;

  // Function that prints the values of the gradient of the basis functions over
  // the face.
  virtual void printBasisDer(std::ostream& out) const = 0;

  // Default virtual destructor.
  virtual ~FeFace() = default;

protected:
  const FaceAbs& face_;
  unsigned dof_;
  const std::vector<std::array<unsigned, 3>>& basisComposition_;
  const QuadRule2D& triaRule_;

  std::vector<Real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

// Penalty parameter equal to max(order^2 / diameter_k) for k = k1, k2,
// with order the order of the polynomials k1 and k2 the two polyhedra sharing
// the face, diameter_k their diameter. For external face there is only one
// polyhedron and so no need for the maximum.
  Real penaltyParam_;

  // Auxiliary function that computes the values of the basis functions and
  // their gradient and fills phi_ and phDer_.
  virtual void compute_basis() = 0;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeFace::getDof() const
{
  return dof_;
}

inline unsigned FeFace::getElemIn() const
{
  return face_.getTetIn().getPoly().getId();
}

inline SizeType FeFace::getQuadPointsNo() const
{
  return triaRule_.getPointsNo();
}

inline Real FeFace::getWeight(SizeType p) const
{
  return triaRule_.getWeight(p);
}

inline Eigen::Vector3d FeFace::getQuadPoint(SizeType p) const
{
  return face_.getTetIn().getMap() * (QuadRuleManager::instance().getFaceMap(face_.getFaceNoTetIn()) *
                                     triaRule_.getPoint(p).homogeneous());
}

inline Real FeFace::getPenaltyParam() const
{
  return penaltyParam_;
}

inline Real FeFace::getAreaDoubled() const
{
  return face_.getAreaDoubled();
}

inline const Eigen::Vector3d& FeFace::getNormal() const
{
  return face_.getNormal();
}

} // namespace PolyDG

#endif // _FE_FACE_HPP_
