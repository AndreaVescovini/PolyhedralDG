#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "MeshProxy.hpp"
#include "MeshReader.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"
#include "FaceExt.hpp"
#include "FaceInt.hpp"
#include "Polyhedron.hpp"

namespace dgfem {

using geom::Vertex;
using geom::Tetrahedron;
using geom::Polyhedron;
using geom::FaceExt;
using geom::FaceInt;

class MeshReader;

class Mesh
{
public:
  Mesh() = default;

  Mesh(const std::string& fileName, MeshReader& reader);

  void printAll(std::ostream& out = std::cout) const;
  void printHead(std::ostream& out = std::cout) const;

  friend class MeshProxy;

private:
  std::vector<Vertex> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  std::vector<FaceExt> facesExt_;
  std::vector<FaceInt> facesInt_;
  std::vector<Polyhedron> polyhedra_;
  // unsigned verticesNo_;
  // unsigned tetrahedraNo_;
  // unsigned polyhedraNo_;
  // unsigned facesNo_;

  // void computePolyInfo();
  void print(unsigned lineNo, std::ostream& out = std::cout) const;
};

}

#endif // _MESH_HPP_
