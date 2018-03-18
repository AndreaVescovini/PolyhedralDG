#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FeFace.hpp"
#include "FaceExt.hpp"
#include "QuadRuleManager.hpp"
#include <vector>
#include <array>
#include "geom.hpp"
#include <Eigen/Dense>

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

  void printBasis(std::ostream& out) const override;
  void printBasisDer(std::ostream& out) const override;

  virtual ~FeFaceExt() = default;

private:
  const TheFace& face_;

  void compute_basis() override;
  unsigned sub2ind(unsigned p, unsigned f) const;

};

}

#endif // _FE_FACE_EXT_HPP_
