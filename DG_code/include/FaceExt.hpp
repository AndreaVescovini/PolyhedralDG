#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include <array>
#include <iostream>
#include "Face.hpp"

namespace geom
{

class FaceExt : public Face
{
public:
  // FaceExt() = default;
  FaceExt(const Vertex& v1, const Vertex& v2, const Vertex& v3, unsigned BClabel);

  unsigned getBClabel() const;

private:
  unsigned BClabel_; // mettere l'enum

  void print(std::ostream& out) const override;
};

}

#endif // _FACE_EXT_HPP_
