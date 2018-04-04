#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FeFace.hpp"
#include "FaceExt.hpp"
#include "QuadRuleManager.hpp"
#include <vector>
#include <array>
#include "geom.hpp"
#include <Eigen/Core>

namespace dgfem
{

class FeFaceExt : public FeFace
{
public:
  using TheFace = geom::FaceExt;

  FeFaceExt(const TheFace& face, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRuleManager::Rule2D& triaRule);

  geom::real getPhi(unsigned p, unsigned f) const;
  const Eigen::Vector3d& getPhiDer(unsigned p, unsigned f) const;

  unsigned getElem() const;
  unsigned getBClabel() const;

  Eigen::Vector3d getQuadPoint(unsigned q) const override;

  geom::real getAreaDoubled() const override;
  const Eigen::Vector3d& getNormal() const override;

  void printBasis(std::ostream& out) const override;
  void printBasisDer(std::ostream& out) const override;

  virtual ~FeFaceExt() = default;

private:
  const TheFace& face_;

// Auxiliary function that fills phi_ and phiDer_
  void compute_basis() override;

// Auxiliary function that, given the quadrature point p and basis function f,
// returns the index in which the corresponding value is stored in phi_ and phiDer_
  unsigned sub2ind(unsigned p, unsigned f) const;

};

} // namespace dgfem

#endif // _FE_FACE_EXT_HPP_
