#ifndef _FACE_HPP_
#define _FACE_HPP_

#include <array>
#include <iostream>
#include "geom.hpp"

namespace geom
{

class Face
{
public:
  Face() = default;
  explicit Face(const std::array<labelType, 3> vertices);

  labelType getTet1() const;
  labelType getFtet1() const;

  real getArea() const;
  // get normal vector

  friend std::ostream& operator<<(std::ostream& out, const Face& face);

protected:
  std::array<labelType, 3> vertices_;
  labelType tet1_;
  labelType ftet1_;
  real area_;
  // normal vector

  virtual void print(std::ostream& out) const = 0;

};

}

#endif // _FACE_HPP_
