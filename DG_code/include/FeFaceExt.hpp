#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FaceExt.hpp"

namespace dgfem
{

class FeFaceExt
{
public:
  using theFaceExt = geom::FaceExt;

  explicit FeFaceExt(const theFaceExt& face);

  virtual ~FeFaceExt() = default;
  
private:
  const theFaceExt& face_;
};

}

#endif // _FE_FACE_EXT_HPP_
