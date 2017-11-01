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

class MeshReader;

class Mesh
{
public:
  Mesh() = default;

  Mesh(const std::string& fileName, MeshReader& reader);

  // void load(const std::string& fileName);

  void print(unsigned lineNo, std::ostream& out = std::cout) const;
  void printAll(std::ostream& out) const;
  void printHead(std::ostream& out = std::cout) const;

  // void setVertices(const std::vector<Point>& vertices);
  // void setTetrahedra(const std::vector<Tetrahedron>& tetrahedra);
  // void setPolyhedra(const std::vector<Polyhedron>& polyhedra);
  // void setVerticesNo(unsigned verticesNo);
  // void setTetrahedraNo(unsigned tetrahedraNo);
  // void setPolyhedraNo(unsigned polyhedraNo);

  friend class MeshProxy;

private:
  std::vector<geom::Point> vertices_;
  std::vector<geom::Tetrahedron> tetrahedra_;
  std::vector<geom::Polyhedron> polyhedra_;
  std::vector<geom::Face> faces_;
  // unsigned verticesNo_;
  // unsigned tetrahedraNo_;
  // unsigned polyhedraNo_;
  // unsigned facesNo_;
};

}

#endif // _MESH_HPP_
