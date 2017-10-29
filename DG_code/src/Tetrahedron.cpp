#include "Tetrahedron.hpp"

namespace geom {

Tetrahedron::Tetrahedron(const std::array<labelType, 4>& vertices, labelType polyNo)
  : vertices_{vertices}, polyNo_{polyNo} {}

std::array<labelType, 4> Tetrahedron::getVertices() const
{
  return vertices_;
}

labelType Tetrahedron::getPolyNo() const
{
  return polyNo_;
}

void Tetrahedron::setPolyNo(labelType polyNo)
{
  polyNo_ = polyNo;
}

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra)
{
  out << "V: " << tetra.vertices_[0] << " " << tetra.vertices_[1] << " "
               << tetra.vertices_[2] << " " << tetra.vertices_[3] << ", P: "
               << tetra.polyNo_;
  return out;
}

}
