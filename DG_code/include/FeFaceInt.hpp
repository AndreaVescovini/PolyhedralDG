#ifndef _FE_FACE_INT_HPP_
#define _FE_FACE_INT_HPP_

#include "FaceInt.hpp"

namespace dgfem
{

class FeFaceInt
{
public:
  using theFaceInt = geom::FaceInt;

  explicit FeFaceInt(const theFaceInt& face);

  virtual ~FeFaceInt() = default;
private:
  const theFaceInt& face_;
};

}

#endif // _FE_FACE_INT_HPP_
