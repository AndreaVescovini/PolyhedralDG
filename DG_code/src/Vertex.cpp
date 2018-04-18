#include "Vertex.hpp"

// #include <cmath>

namespace PolyDG
{

Vertex::Vertex(Real x, Real y, Real z)
  : id_{counter_}, coords_{x, y, z}
{
  counter_++;
}

// void Vertex::setCoords(const Eigen::Vector3d& coords)
// {
//   coords_ = coords;
// }

// Real Vertex::distance(const Vertex& v2) const
// {
//   return std::sqrt( (this->coords_[0] - v2.coords_[0]) * (this->coords_[0] - v2.coords_[0])
//                   + (this->coords_[1] - v2.coords_[1]) * (this->coords_[1] - v2.coords_[1])
//                   + (this->coords_[2] - v2.coords_[2]) * (this->coords_[2] - v2.coords_[2]) );
// }

std::ostream& operator<<(std::ostream& out, const Vertex& v)
{
  out << v.id_ << " " << v.coords_[0] << " " << v.coords_[1] << " " << v.coords_[2];
  return out;
}

// bool compX(const Vertex& lhs, const Vertex& rhs)
// {
//   return lhs.coords_[0] < rhs.coords_[0];
// }
//
// bool compY(const Vertex& lhs, const Vertex& rhs)
// {
//   return lhs.coords_[1] < rhs.coords_[1];
// }
//
// bool compZ(const Vertex& lhs, const Vertex& rhs)
// {
//   return lhs.coords_[2] < rhs.coords_[2];
// }

unsigned Vertex::counter_ = 0;

} // namespace PolyDG
