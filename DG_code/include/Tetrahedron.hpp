#ifndef _TETRAHEDRON_HPP_
#define _TETRAHEDRON_HPP_

#include <array>
#include <iostream>

namespace geom {

using labelType = unsigned; // da includere poi da qualche altra parte

class Tetrahedron
{
public:
  Tetrahedron() = default;
  Tetrahedron(const std::array<labelType, 4>& vertices, labelType polyNo = 0);

  std::array<labelType, 4> getVertices() const;
  labelType getPolyNo() const;
  void setPolyNo(labelType polyNo);

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  std::array<labelType, 4> vertices_;
  labelType polyNo_;
};

}

#endif // _TETRAHEDRON_HPP_
