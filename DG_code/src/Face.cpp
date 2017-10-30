#include "Face.hpp"

namespace geom {

labelType Face::getTet1() const
{
  return tet1_;
}

labelType Face::getTet2() const
{
  return tet2_;
}

labelType Face::getFtet1() const
{
  return ftet1_;
}

labelType Face::getFtet2() const
{
  return ftet2_;
}

real Face::getArea() const
{
  return area_;
}

}
