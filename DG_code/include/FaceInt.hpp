#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include "Face.hpp"
#include "geom.hpp"

namespace geom
{

class FaceInt : public Face
{
public:
  FaceInt() = default;

  labelType getTet2() const;
  labelType getFtet2() const;

private:
  labelType tet2_;
  labelType ftet2_;

};

}

#endif // _FACE_INT_HPP_
