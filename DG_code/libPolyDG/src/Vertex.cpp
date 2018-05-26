/*!
    @file   Vertex.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class Vertex
*/

#include "Vertex.hpp"

namespace PolyDG
{

Vertex::Vertex(Real x, Real y, Real z)
  : id_{counter_}, coords_{x, y, z}
{
  counter_++;
}

std::ostream& operator<<(std::ostream& out, const Vertex& v)
{
  out << v.id_ << " " << v.coords_(0) << " " << v.coords_(1) << " " << v.coords_(2);
  return out;
}

unsigned Vertex::counter_ = 0;

} // namespace PolyDG
