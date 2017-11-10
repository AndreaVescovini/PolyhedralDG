#include "Polyhedron.hpp"

namespace geom {

Polyhedron::Polyhedron()
  : id_{counter_}
{
  counter_++;
}

real Polyhedron::getDiameter() const
{
  return diameter_;
}

void Polyhedron::addTetra(const Tetrahedron& tet)
{
  tetrahedra_.emplace_back(tet);
}

const Tetrahedron& Polyhedron::getTetra(unsigned i) const
{
  return tetrahedra_[i];
}

unsigned Polyhedron::getId() const
{
  return id_;
}

// void Polyhedron::computeBBandDiam()
// {
//   for(const Tetrahedron& t : tetrahedra_)
//     t.getVerticesCoords();
// }

void Polyhedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

unsigned Polyhedron::counter_ = 0;

}
