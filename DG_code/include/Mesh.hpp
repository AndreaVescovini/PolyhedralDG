#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "MeshProxy.hpp"
#include "MeshReader.hpp"
#include "Point.hpp"
#include "Tetrahedron.hpp"
#include "Polyhedron.hpp"
#include "Face.hpp"

namespace dgfem {

using geom::Point;
using geom::Tetrahedron;
using geom::Polyhedron;
using geom::Face;

class MeshReader;

class Mesh
{
public:
  Mesh() = default;

  Mesh(const std::string& fileName, MeshReader& reader);

  void print(unsigned lineNo, std::ostream& out = std::cout) const;
  void printAll(std::ostream& out) const;
  void printHead(std::ostream& out = std::cout) const;

  friend class MeshProxy;

private:
  std::vector<Point> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  std::vector<Polyhedron> polyhedra_;
  std::vector<Face> faces_;
  // unsigned verticesNo_;
  // unsigned tetrahedraNo_;
  // unsigned polyhedraNo_;
  // unsigned facesNo_;
};

}

#endif // _MESH_HPP_
