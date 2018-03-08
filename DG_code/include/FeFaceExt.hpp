#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FaceExt.hpp"
#include "QuadRule.hpp"

namespace dgfem
{

class FeFaceExt
{
public:
  using theFaceExt = geom::FaceExt;

  explicit FeFaceExt(const theFaceExt& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule);

  unsigned getOrder() const;
  unsigned getDofNo() const;

  virtual ~FeFaceExt() = default;
private:
  const theFaceExt& face_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
};

}

#endif // _FE_FACE_EXT_HPP_
