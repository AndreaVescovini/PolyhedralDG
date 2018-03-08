#ifndef _FE_FACE_INT_HPP_
#define _FE_FACE_INT_HPP_

#include "FaceInt.hpp"
#include "QuadRule.hpp"

namespace dgfem
{

class FeFaceInt
{
public:
  using theFaceInt = geom::FaceInt;

  explicit FeFaceInt(const theFaceInt& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule);

  unsigned getOrder() const;
  unsigned getDofNo() const;

  virtual ~FeFaceInt() = default;
private:
  const theFaceInt& face_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
};

}

#endif // _FE_FACE_INT_HPP_
