#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "Point.hpp"
#include "Tetrahedron.hpp"
#include "Polyhedron.hpp"
#include "Face.hpp"

namespace geom {

class Mesh
{
public:
  Mesh() = default;

  void load(const std::string& fileName);
  
  void print(unsigned lineNo, std::ostream& out = std::cout) const;
  void printAll(std::ostream& out) const;
  void printHead(std::ostream& out = std::cout) const;

  void setVertices(const std::vector<Point>& vertices);
  void setTetrahedra(const std::vector<Tetrahedron>& tetrahedra);
  void setPolyhedra(const std::vector<Polyhedron>& polyhedra);
  void setVerticesNo(unsigned verticesNo);
  void setTetrahedraNo(unsigned tetrahedraNo);
  void setPolyhedraNo(unsigned polyhedraNo);

private:
  std::vector<Point> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  std::vector<Polyhedron> polyhedra_;
  std::vector<Face> faces_;
  unsigned verticesNo_;
  unsigned tetrahedraNo_;
  unsigned polyhedraNo_;
  unsigned facesNo_;
};

}

#endif // _MESH_HPP_
