#ifndef _FACE_HPP_
#define _FACE_HPP_

#include <array>
#include "geom.hpp"

namespace geom
{

class Face
{
public:
  Face() = default;

  labelType getTet1() const;
  labelType getTet2() const;
  labelType getFtet1() const;
  labelType getFtet2() const;

  real getArea() const;


private:
  std::array<labelType, 3> vertices_;
  labelType tet1_;
  labelType tet2_;
  labelType ftet1_;
  labelType ftet2_;
  real area_;
  // normal vector

};

}

#endif // _FACE_HPP_
