/*!
    @file   Polyhedron.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class Polyhedron
*/

#include "Polyhedron.hpp"

#include <algorithm>
#include <iterator>

namespace PolyDG
{

Polyhedron::Polyhedron()
  : id_{counter_}, diameter_{0.0}
{
  counter_++;
}

void Polyhedron::computeBB()
{
  // I store in boundingBox_ the min and max vertices coordinates.
  for(auto it = verticesExt_.cbegin(); it != verticesExt_.cend(); it++)
    boundingBox_.extend(it->get().getCoords());
}

void Polyhedron::computeDiameter()
{
  diameter_ = 0.0;

  for(auto i = verticesExt_.cbegin(); i != verticesExt_.cend(); i++)
    for(auto j = std::next(i); j != verticesExt_.cend(); j++)
      diameter_ = std::max(diameter_, i->get().distance(*j));
}

std::ostream& operator<<(std::ostream& out, const Polyhedron& poly)
{
  out << poly.id_ << " " << "T:";
  for(const Tetrahedron& tet : poly.tetrahedra_)
    out << ' ' << tet.getId();

  return out;
}

unsigned Polyhedron::counter_ = 0;

} // namespace PolyDG
