#ifndef _FE_FACE_HPP_
#define _FE_FACE_HPP_

#include "PolyDG.hpp"
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
  FeFace(unsigned dof, const std::vector<std::array<unsigned, 3>>& basisComposition,
         const QuadRuleManager::Rule2D& triaRule);

  inline unsigned getDof() const;

  inline SizeType getQuadPointsNo() const;
  inline Real getWeight(SizeType p) const;
  virtual Eigen::Vector3d getQuadPoint(SizeType p) const = 0;

  inline Real getPenaltyParam() const;

  virtual Real getAreaDoubled() const = 0;
  virtual const Eigen::Vector3d& getNormal() const = 0;

  virtual void printBasis(std::ostream& out) const = 0;
  virtual void printBasisDer(std::ostream& out) const = 0;

  virtual ~FeFace() = default;

protected:
  unsigned dof_;
  const std::vector<std::array<unsigned, 3>>& basisComposition_;
  const QuadRuleManager::Rule2D& triaRule_;

  std::vector<Real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

// Penalty parameter equal to max(order^2 / diameter_k) for k = k1, k2,
// with order the order of the polynomials k1 and k2 the two polyhedra sharing
// the face, diameter_k their diameter. For external face there is only one
// polyhedron and so no need for the maximum.
  Real penaltyParam_;

  virtual void compute_basis() = 0;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeFace::getDof() const
{
  return dof_;
}

inline SizeType FeFace::getQuadPointsNo() const
{
  return triaRule_.getPointsNo();
}

inline Real FeFace::getWeight(SizeType p) const
{
  return triaRule_.getWeight(p);
}

inline Real FeFace::getPenaltyParam() const
{
  return penaltyParam_;
}

} // namespace PolyDG

#endif // _FE_FACE_HPP_
