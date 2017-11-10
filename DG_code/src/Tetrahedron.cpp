#include "Tetrahedron.hpp"

namespace geom {

Tetrahedron::Tetrahedron(const Vertex& v1, const Vertex& v2, const Vertex& v3,
                         const Vertex& v4, const Polyhedron* poly)
  : id_{counter_}, vertices_{{v1, v2, v3, v4}}, poly_{poly}
{
  counter_++;
}

Tetrahedron::Tetrahedron(const Vertex& v1, const Vertex& v2, const Vertex& v3,
                         const Vertex& v4, const Polyhedron& poly)
  : id_{counter_}, vertices_{{v1, v2, v3, v4}}, poly_{&poly}
{
  counter_++;
}

const Vertex& Tetrahedron::getVertex(unsigned i) const
{
  return vertices_[i];
}

const Polyhedron& Tetrahedron::getPoly() const
{
  return *poly_;
}

void Tetrahedron::setPoly(const Polyhedron* poly)
{
  poly_ = poly;
}

void Tetrahedron::setPoly(const Polyhedron& poly)
{
  poly_ = &poly;
}

unsigned Tetrahedron::getId() const
{
  return id_;
}

void Tetrahedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra)
{
  out << tetra.id_ << " " << "V: " << tetra.vertices_[0].get().getId() << " "
                                   << tetra.vertices_[1].get().getId() << " "
                                   << tetra.vertices_[2].get().getId() << " "
                                   << tetra.vertices_[3].get().getId() << ", P: "
                                   << tetra.poly_->getId();
  return out;
}

unsigned Tetrahedron::counter_ = 0;

}
