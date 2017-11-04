#include "Face.hpp"

namespace geom {

Face::Face(const std::array<labelType, 3> vertices)
  : vertices_{vertices} {}

labelType Face::getTet1() const
{
  return tet1_;
}

labelType Face::getFtet1() const
{
  return ftet1_;
}

real Face::getArea() const
{
  return area_;
}

std::ostream& operator<<(std::ostream& out, const Face& face)
{
  face.print(out);
  return out;
}

}
