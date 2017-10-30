#include "Tetrahedron.hpp"

namespace geom {

Tetrahedron::Tetrahedron(const std::array<labelType, 4>& vertices, labelType poly)
  : vertices_{vertices}, poly_{poly} {}

std::array<labelType, 4> Tetrahedron::getVertices() const
{
  return vertices_;
}

labelType Tetrahedron::getPoly() const
{
  return poly_;
}

void Tetrahedron::setPoly(labelType poly)
{
  poly_ = poly;
}

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra)
{
  out << "V: " << tetra.vertices_[0] << " " << tetra.vertices_[1] << " "
               << tetra.vertices_[2] << " " << tetra.vertices_[3] << ", P: "
               << tetra.poly_;
  return out;
}

}
