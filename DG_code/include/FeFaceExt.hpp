#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FaceExt.hpp"
#include "FeFace.hpp"
#include "PolyDG.hpp"
#include "QuadRuleManager.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace PolyDG
{

class FeFaceExt : public FeFace
{
public:
  using TheFace = FaceExt;

  FeFaceExt(const TheFace& face, unsigned order, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule2D& triaRule);

  inline Real getPhi(SizeType p, SizeType f) const;
  inline const Eigen::Vector3d& getPhiDer(SizeType p, SizeType f) const;

  inline unsigned getElem() const;
  inline BCType getBClabel() const;

  inline Eigen::Vector3d getQuadPoint(SizeType p) const override;

  inline Real getAreaDoubled() const override;
  inline const Eigen::Vector3d& getNormal() const override;

  void printBasis(std::ostream& out) const override;
  void printBasisDer(std::ostream& out) const override;

  virtual ~FeFaceExt() = default;

private:
  const TheFace& face_;

// Auxiliary function that fills phi_ and phiDer_
  void compute_basis() override;

// Auxiliary function that, given the quadrature point p and basis function f,
// returns the index in which the corresponding value is stored in phi_ and phiDer_
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

inline unsigned FeFaceExt::getElem() const
{
  return face_.getTetIn().getPoly().getId();
}

inline BCType FeFaceExt::getBClabel() const
{
  return face_.getBClabel();
}

inline Eigen::Vector3d FeFaceExt::getQuadPoint(SizeType p) const
{
  return face_.getTetIn().getMap() * (QuadRuleManager::getFaceMap(face_.getFaceNoTetIn()) *
                                     triaRule_.getPoint(p).homogeneous());
}

inline Real FeFaceExt::getAreaDoubled() const
{
  return face_.getAreaDoubled();
}

inline const Eigen::Vector3d& FeFaceExt::getNormal() const
{
  return face_.getNormal();
}

inline SizeType FeFaceExt::sub2ind(SizeType p, SizeType f) const
{
  return f + dof_ * p;
}

} // namespace PolyDG

#endif // _FE_FACE_EXT_HPP_
