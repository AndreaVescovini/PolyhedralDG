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

void Polyhedron::computeBB()
{
  boundingBox_[0][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compX)->get().getX();
  boundingBox_[0][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compX)->get().getX();
  boundingBox_[1][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compY)->get().getY();
  boundingBox_[1][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compY)->get().getY();
  boundingBox_[2][0] = std::min_element(vertices_.cbegin(), vertices_.cend(), compZ)->get().getZ();
  boundingBox_[2][1] = std::max_element(vertices_.cbegin(), vertices_.cend(), compZ)->get().getZ();
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
