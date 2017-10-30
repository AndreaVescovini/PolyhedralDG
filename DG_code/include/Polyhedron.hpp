#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include <array>
#include "geom.hpp"

namespace geom
{

class Polyhedron
{
public:
  Polyhedron() = default;

  real getDiameter() const;

private:
  real diameter_;
  std::array<std::array<real, 2>, 3> boundingBox_;
};

}

#endif // _POLYHEDRON_HPP_
