#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include <vector>
#include <array>
#include <functional>
#include <unordered_set>
#include "geom.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"

namespace geom
{

class Tetrahedron;

class Polyhedron
{
public:
  Polyhedron();

  real getDiameter() const;
  void addTetra(const Tetrahedron& tet);
  const Tetrahedron& getTetra(unsigned i) const;
  unsigned getId() const;
  void computeBB();
  void computeDiameter();

  static void resetCounter(unsigned counter = 0);

private:
  const unsigned id_;
  std::vector<std::reference_wrapper<const Tetrahedron>> tetrahedra_;
  // std::unordered_set<std::reference_wrapper<const Vertex>> vertices_;
  std::unordered_set<std::reference_wrapper<const Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> vertices_;
  std::array<interval, 3> boundingBox_;
  real diameter_;

  static unsigned counter_;
};

}

#endif // _POLYHEDRON_HPP_
