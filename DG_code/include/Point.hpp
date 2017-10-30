#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <array>
#include <iostream>
#include "geom.hpp"

namespace geom {

class Point
{
public:
  Point() = default;
  explicit Point(const std::array<real, 3>& coords);

  std::array<real, 3> getCoords() const;
  real getX() const;
  real getY() const;
  real getZ() const;

  real distance(const Point& p2) const;

  friend std::ostream& operator<<(std::ostream& out, const Point& point);

private:
  std::array<real, 3> coords_;
};

}

#endif // _POINT_HPP_
