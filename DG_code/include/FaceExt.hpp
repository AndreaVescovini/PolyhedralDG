#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include <array>
#include <iostream>
#include "Face.hpp"
#include "geom.hpp"

namespace geom
{

class FaceExt : public Face
{
public:
  FaceExt() = default;
  FaceExt(const std::array<labelType, 3> vertices, unsigned BClabel);

  unsigned getBClabel() const;

private:
  unsigned BClabel_;

  void print(std::ostream& out) const;

};

}

#endif // _FACE_EXT_HPP_
