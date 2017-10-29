#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>
#include <string>
#include <iostream>

#include "Point.hpp"
#include "Tetrahedron.hpp"
// #include "MeshReader.hpp"

namespace geom {

class Mesh
{
public:
  Mesh() = default;

  void load(const std::string& fileName);
  void print(unsigned lineNo, std::ostream& out = std::cout) const;
  void printAll(std::ostream& out) const;
  void printHead(std::ostream& out = std::cout) const;

private:
  std::vector<Point> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  unsigned verticesNo_;
  unsigned tetrahedraNo_;
  unsigned polyhedraNo_;
};

}

#endif // _MESH_HPP_
