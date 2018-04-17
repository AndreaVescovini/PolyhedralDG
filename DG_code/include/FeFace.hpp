#ifndef _FE_FACE_
#define _FE_FACE_

#include "QuadRuleManager.hpp"
#include <vector>
#include <array>
#include "geom.hpp"
#include <Eigen/Core>
#include <iostream>

namespace dgfem
{

class FeFace
{
public:
  FeFace(unsigned dofNo, const std::vector<std::array<unsigned, 3>>& basisComposition,
         const QuadRuleManager::Rule2D& triaRule);

  inline unsigned getDofNo() const;

  inline unsigned getQuadPointsNo() const;
  inline geom::real getWeight(unsigned q) const;
  virtual Eigen::Vector3d getQuadPoint(unsigned q) const = 0;

  inline geom::real getPenaltyParam() const;

  virtual geom::real getAreaDoubled() const = 0;
  virtual const Eigen::Vector3d& getNormal() const = 0;

  virtual void printBasis(std::ostream& out) const = 0;
  virtual void printBasisDer(std::ostream& out) const = 0;

  virtual ~FeFace() = default;

protected:
  unsigned dofNo_;
  const std::vector<std::array<unsigned, 3>>& basisComposition_;
  const QuadRuleManager::Rule2D& triaRule_;

  std::vector<geom::real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

// Penalty parameter equal to max(order^2 / diameter_k) for k = k1, k2,
// with order the order of the polynomials k1 and k2 the two polyhedra sharing
// the face, diameter_k their diameter. For external face there is only one
// polyhedron and so no need for the maximum.
  geom::real penaltyParam_;

  virtual void compute_basis() = 0;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeFace::getDofNo() const
{
  return dofNo_;
}

inline unsigned FeFace::getQuadPointsNo() const
{
  return triaRule_.getPointsNo();
}

inline geom::real FeFace::getWeight(unsigned q) const
{
  return triaRule_.getWeight(q);
}

inline geom::real FeFace::getPenaltyParam() const
{
  return penaltyParam_;
}

} // namespace dgfem

#endif // _FE_FACE_
