#ifndef _FE_FACE_INT_HPP_
#define _FE_FACE_INT_HPP_

#include "FaceInt.hpp"
#include "FeFace.hpp"
#include "PolyDG.hpp"
#include "QuadRuleManager.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace PolyDG
{

class FeFaceInt : public FeFace
{
public:
  using TheFace = FaceInt;

  FeFaceInt(const TheFace& face, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule2D& triaRule);

  inline Real getPhi(SideType s, unsigned p, unsigned f) const;
  inline const Eigen::Vector3d& getPhiDer(SideType s, unsigned p, unsigned f) const;

  inline unsigned getElem(SideType s) const;

  inline Eigen::Vector3d getQuadPoint(unsigned q) const override;

  inline Real getAreaDoubled() const override;
  inline const Eigen::Vector3d& getNormal() const override;

  void printBasis(std::ostream& out) const override;
  void printBasisDer(std::ostream& out) const override;

  virtual ~FeFaceInt() = default;

private:
  const TheFace& face_;

// Auxiliary function that fills phi_ and phiDer_
  void compute_basis() override;

// Auxiliary function that, given the s of the face, the quadrature point p
// and basis function f, returns the index in which the corresponding value
// is stored in phi_ and phiDer_
  inline unsigned sub2ind(SideType s, unsigned p, unsigned f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceInt::getPhi(SideType s, unsigned p, unsigned f) const
{
  return phi_[sub2ind(s, p, f)];
}

inline const Eigen::Vector3d& FeFaceInt::getPhiDer(SideType s, unsigned p, unsigned f) const
{
  return phiDer_[sub2ind(s, p, f)];
}

inline unsigned FeFaceInt::getElem(SideType s) const
{
  return (s == Out ? face_.getTetIn().getPoly().getId() : face_.getTetOut().getPoly().getId());
}

inline Eigen::Vector3d FeFaceInt::getQuadPoint(unsigned q) const
{
  return face_.getTetIn().getMap() * (QuadRuleManager::getFaceMap(face_.getFaceNoTetIn()) *
                                     triaRule_.getPoint(q).homogeneous());
}

inline Real FeFaceInt::getAreaDoubled() const
{
  return face_.getAreaDoubled();
}

inline const Eigen::Vector3d& FeFaceInt::getNormal() const
{
  return face_.getNormal();
}

inline unsigned FeFaceInt::sub2ind(SideType s, unsigned p, unsigned f) const
{
  return (s == In) + 2 * (f + p * dofNo_);
}

} // namespace PolyDG

#endif // _FE_FACE_INT_HPP_
