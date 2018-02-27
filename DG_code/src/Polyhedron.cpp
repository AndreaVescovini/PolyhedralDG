#include "Polyhedron.hpp"
#include <algorithm>
#include <iterator>

namespace geom {

Polyhedron::Polyhedron()
  : id_{counter_}, diameter_{0.0}
{
  counter_++;
}

real Polyhedron::getDiameter() const
{
  return diameter_;
}

void Polyhedron::addTetra(Tetrahedron& tet)
{
  tetrahedra_.emplace_back(tet);
}

void Polyhedron::addVertex(Vertex& v)
{
  vertices_.emplace(v);
}

const Tetrahedron& Polyhedron::getTetra(unsigned i) const
{
  return tetrahedra_[i];
}

unsigned Polyhedron::getId() const
{
  return id_;
}

void Polyhedron::computeBB()
{
  // I store in boundingBox_ the min and max vertices coordinates.
  // boundingBox_[0][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compX)->get().getX();
  // boundingBox_[0][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compX)->get().getX();
  // boundingBox_[1][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compY)->get().getY();
  // boundingBox_[1][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compY)->get().getY();
  // boundingBox_[2][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compZ)->get().getZ();
  // boundingBox_[2][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compZ)->get().getZ();

  for(auto it = vertices_.cbegin(); it != vertices_.cend(); it++)
    boundingBox_.extend(it->get().getCoords());
}

void Polyhedron::computeDiameter()
{
  diameter_ = 0.0;

  for(auto i = vertices_.cbegin(); i != vertices_.cend(); i++)
    for(auto j = std::next(i); j != vertices_.cend(); j++)
      diameter_ = std::max(diameter_, i->get().distance(*j));
}

void Polyhedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

unsigned Polyhedron::counter_ = 0;

}
