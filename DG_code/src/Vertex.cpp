#include "Vertex.hpp"
#include <cmath>

namespace geom {

// Vertex::Vertex()
//   : id_{counter_}, coords_{{0.0,0.0,0.0}}
// {
//   counter_++;
// }

Vertex::Vertex(real x, real y, real z)
  : id_{counter_}, coords_{x, y, z}
{
  counter_++;
}

const Eigen::Vector3d& Vertex::getCoords() const
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

// real Vertex::distance(const Vertex& v2) const
// {
//   return std::sqrt( (this->coords_[0] - v2.coords_[0]) * (this->coords_[0] - v2.coords_[0])
//                   + (this->coords_[1] - v2.coords_[1]) * (this->coords_[1] - v2.coords_[1])
//                   + (this->coords_[2] - v2.coords_[2]) * (this->coords_[2] - v2.coords_[2]) );
// }

real Vertex::distance(const Vertex& v2) const
{
  return (this->getCoords() - v2.getCoords()).norm();
}

std::ostream& operator<<(std::ostream& out, const Vertex& v)
{
  out << v.id_ << " " << v.coords_[0] << " " << v.coords_[1] << " " << v.coords_[2];
  return out;
}

bool compX(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.coords_[0] < rhs.coords_[0];
}

bool compY(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.coords_[1] < rhs.coords_[1];
}

bool compZ(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.coords_[2] < rhs.coords_[2];
}

void Vertex::resetCounter(unsigned counter)
{
  counter_ = counter;
}

unsigned Vertex::counter_ = 0;

}
