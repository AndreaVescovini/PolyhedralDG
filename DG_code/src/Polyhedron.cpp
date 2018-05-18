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

void Polyhedron::addTetra(Tetrahedron& tet)
{
  tetrahedra_.emplace_back(tet);

  for(SizeType i = 0; i < 4; i++)
  {
    auto insertion = vertices_.emplace(tet.getVertex(i));
    if(insertion.second == true)
    {
      boundingBox_.extend(tet.getVertex(i).getCoords());
      for(auto iter = vertices_.cbegin(); iter != vertices_.cend(); iter++)
        diameter_ = std::max(diameter_, iter->get().distance(tet.getVertex(i)));
    }
  }
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
