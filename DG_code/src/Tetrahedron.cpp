#include "Tetrahedron.hpp"

namespace geom {

Tetrahedron::Tetrahedron(Vertex& v1,  Vertex& v2, Vertex& v3, Vertex& v4,
                         Polyhedron* poly)
  : id_{counter_}, vertices_{{v1, v2, v3, v4}}, poly_{poly}
{
  counter_++;
}

Tetrahedron::Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4,
                         Polyhedron& poly)
  : Tetrahedron(v1, v2, v3, v4, &poly) {}

const Vertex& Tetrahedron::getVertex(unsigned i) const
{
  return vertices_[i];
}

Vertex& Tetrahedron::getVertex(unsigned i)
{
  return vertices_[i];
}

const Polyhedron& Tetrahedron::getPoly() const
{
  return *poly_;
}

Polyhedron& Tetrahedron::getPoly()
{
  return *poly_;
}

void Tetrahedron::setPoly(Polyhedron* poly)
{
  poly_ = poly;
}

void Tetrahedron::setPoly(Polyhedron& poly)
{
  setPoly(&poly);
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
