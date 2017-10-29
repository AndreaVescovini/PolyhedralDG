#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <array>
#include <iostream>

namespace geom {

using real = double; // da includere poi da qualche altra parte

class Point
{
public:
  Point() = default;
  explicit Point(const std::array<real, 3>& coords);

  std::array<real, 3> getCoords() const;
  real getX() const;
  real getY() const;
  real getZ() const;

  friend std::ostream& operator<<(std::ostream& out, const Point& point);

private:
  std::array<real, 3> coords_;
};



}

#endif // _POINT_HPP_
