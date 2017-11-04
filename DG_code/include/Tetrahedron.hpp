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

  Tetrahedron(const Tetrahedron&) = default;
  Tetrahedron& operator=(const Tetrahedron&) = default;
  Tetrahedron(Tetrahedron&&) = default;
  Tetrahedron& operator=(Tetrahedron&&) = default;

  std::array<labelType, 4> getVertices() const;
  labelType getPoly() const;
  void setPoly(labelType poly);

  virtual ~Tetrahedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  const std::array<labelType, 4> vertices_;
  labelType poly_;
};

}

#endif // _TETRAHEDRON_HPP_
