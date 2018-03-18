#ifndef _FE_FACE_
#define _FE_FACE_

#include "QuadRuleManager.hpp"
#include <vector>
#include <array>
#include "geom.hpp"
#include <Eigen/Dense>
#include <iostream>

namespace dgfem
{

class FeFace
{
public:
  FeFace(unsigned order, unsigned dofNo,
         const std::vector<std::array<unsigned, 3>>& basisComposition,
         const QuadRuleManager::Rule2D& triaRule);

  unsigned getOrder() const;
  unsigned getDofNo() const;

  virtual void printBasis(std::ostream& out) const = 0;
  virtual void printBasisDer(std::ostream& out) const = 0;

  virtual ~FeFace() = default;

protected:
  unsigned order_;
  unsigned dofNo_;
  const std::vector<std::array<unsigned, 3>>& basisComposition_;
  const QuadRuleManager::Rule2D& triaRule_;

  std::vector<geom::real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

  virtual void compute_basis() = 0;

};

}

#endif // _FE_FACE_
