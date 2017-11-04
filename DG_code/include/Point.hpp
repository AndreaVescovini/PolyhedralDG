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

  Point(const Point&) = default;
  Point& operator=(const Point&) = default;
  Point(Point&&) = default;
  Point& operator=(Point&&) = default;

  std::array<real, 3> getCoords() const;
  real getX() const;
  real getY() const;
  real getZ() const;

  real distance(const Point& p2) const;

  virtual ~Point() = default;

  friend std::ostream& operator<<(std::ostream& out, const Point& point);

private:
  const std::array<real, 3> coords_;
};

}

#endif // _POINT_HPP_
