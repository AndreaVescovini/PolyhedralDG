#include "Point.hpp"
#include <cmath>

namespace geom {

Point::Point(const std::array<real, 3>& coords)
  : coords_{coords} {}

std::array<real, 3> Point::getCoords() const
{
  return coords_;
}

real Point::getX() const
{
  return coords_[0];
}

real Point::getY() const
{
  return coords_[1];
}

real Point::getZ() const
{
  return coords_[2];
}

real Point::distance(const Point& p2) const
{
  return std::sqrt( (this->coords_[0] - p2.coords_[0]) * (this->coords_[0] - p2.coords_[0])
                  + (this->coords_[1] - p2.coords_[1]) * (this->coords_[1] - p2.coords_[1])
                  + (this->coords_[2] - p2.coords_[2]) * (this->coords_[2] - p2.coords_[2]) );
}

std::ostream& operator<<(std::ostream& out, const Point& point)
{
  out << point.coords_[0] << " " << point.coords_[1] << " " << point.coords_[2];
  return out;
}

}
