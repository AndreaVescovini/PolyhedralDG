#ifndef _TETRAHEDRON_HPP_
#define _TETRAHEDRON_HPP_

#include <array>
#include <iostream>
#include "geom.hpp"

namespace geom {

class Tetrahedron
{
public:
  Tetrahedron() = default;
  Tetrahedron(const std::array<labelType, 4>& vertices, labelType poly = 0);

  std::array<labelType, 4> getVertices() const;
  labelType getPoly() const;
  void setPoly(labelType poly);

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  std::array<labelType, 4> vertices_;
  labelType poly_;
};

}

#endif // _TETRAHEDRON_HPP_
