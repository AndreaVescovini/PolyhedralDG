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
  // The defolut constructor also sets the id using the counter.
  Polyhedron();

  real getDiameter() const;
  void addTetra(Tetrahedron& tet);
  void addVertex(Vertex& v);
  const Tetrahedron& getTetra(unsigned i) const;
  Tetrahedron& getTetra(unsigned i);
  // getVertex?
  unsigned getId() const;

  // Function that starting from the unordered set of vertices computes the
  // cartesian bounding box of the polyhedron.
  void computeBB();

  // Function that starting from the unordered set of vertices computes the
  // diameter of the polyhderon i.e. the maximum distance between two points
  // inside the polyhderon.
  void computeDiameter();

  static void resetCounter(unsigned counter = 0);

private:
  const unsigned id_;
  std::vector<std::reference_wrapper<Tetrahedron>> tetrahedra_;

  // Here I store vertices coming from faces
  std::unordered_set<std::reference_wrapper<Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> vertices_;
  std::array<interval, 3> boundingBox_;
  real diameter_;

  static unsigned counter_;
};

}

#endif // _POLYHEDRON_HPP_
