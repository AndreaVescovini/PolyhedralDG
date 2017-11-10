#include "Vertex.hpp"
#include <cmath>

namespace geom {

// Vertex::Vertex()
//   : id_{counter_}, coords_{{0.0,0.0,0.0}}
// {
//   counter_++;
// }

Vertex::Vertex(const std::array<real, 3>& coords)
  : id_{counter_}, coords_{coords}
{
  counter_++;
}

std::array<real, 3> Vertex::getCoords() const
{
  return coords_;
}

real Vertex::getX() const
{
  return coords_[0];
}

real Vertex::getY() const
{
  return coords_[1];
}

real Vertex::getZ() const
{
  return coords_[2];
}

unsigned Vertex::getId() const
{
  return id_;
}

real Vertex::distance(const Vertex& v2) const
{
  return std::sqrt( (this->coords_[0] - v2.coords_[0]) * (this->coords_[0] - v2.coords_[0])
                  + (this->coords_[1] - v2.coords_[1]) * (this->coords_[1] - v2.coords_[1])
                  + (this->coords_[2] - v2.coords_[2]) * (this->coords_[2] - v2.coords_[2]) );
}

std::ostream& operator<<(std::ostream& out, const Vertex& v)
{
  out << v.id_ << " " << v.coords_[0] << " " << v.coords_[1] << " " << v.coords_[2];
  return out;
}

void Vertex::resetCounter(unsigned counter)
{
  counter_ = counter;
}

unsigned Vertex::counter_ = 0;

}
