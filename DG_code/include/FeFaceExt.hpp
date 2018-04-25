#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "PolyDG.hpp"
#include "FeFace.hpp"
#include "FaceExt.hpp"
#include "QuadRuleManager.hpp"

#include <Eigen/Core>

#include <vector>
#include <array>

namespace PolyDG
{

class FeFaceExt : public FeFace
{
public:
  using TheFace = FaceExt;

  FeFaceExt(const TheFace& face, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule2D& triaRule);

  inline Real getPhi(unsigned p, unsigned f) const;
  inline const Eigen::Vector3d& getPhiDer(unsigned p, unsigned f) const;

  inline unsigned getElem() const;
  inline BCtype getBClabel() const;

  inline Eigen::Vector3d getQuadPoint(unsigned q) const override;

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
  inline unsigned sub2ind(unsigned p, unsigned f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceExt::getPhi(unsigned p, unsigned f) const
{
  return phi_[sub2ind(p, f)];
}

inline const Eigen::Vector3d& FeFaceExt::getPhiDer(unsigned p, unsigned f) const
{
  return phiDer_[sub2ind(p, f)];
}

inline unsigned FeFaceExt::getElem() const
{
  return face_.getTet1().getPoly().getId();
}

inline BCtype FeFaceExt::getBClabel() const
{
  return face_.getBClabel();
}

inline Eigen::Vector3d FeFaceExt::getQuadPoint(unsigned q) const
{
  return face_.getTet1().getMap() * (QuadRuleManager::getFaceMap(face_.getFaceNoTet1()) *
                                     triaRule_.getPoint(q).homogeneous());
}

inline Real FeFaceExt::getAreaDoubled() const
{
  return face_.getAreaDoubled();
}

inline const Eigen::Vector3d& FeFaceExt::getNormal() const
{
  return face_.getNormal();
}

inline unsigned FeFaceExt::sub2ind(unsigned p, unsigned f) const
{
  return f + dofNo_ * p;
}

} // namespace PolyDG

#endif // _FE_FACE_EXT_HPP_
