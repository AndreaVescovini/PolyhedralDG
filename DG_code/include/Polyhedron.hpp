#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include <vector>
#include <array>
#include <functional>
#include "geom.hpp"
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
  // void computeBBandDiam();

  static void resetCounter(unsigned counter = 0);

private:
  const unsigned id_;
  std::vector<std::reference_wrapper<const Tetrahedron>> tetrahedra_;
  std::array<interval, 3> boundingBox_;
  real diameter_;

  static unsigned counter_;
};

}

#endif // _POLYHEDRON_HPP_
