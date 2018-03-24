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

namespace geom {

class MeshReader;

class Mesh
{
public:
  Mesh() = default;

  Mesh(const std::string& fileName, MeshReader& reader);

  const FaceExt& getFaceExt(unsigned n) const;
  const FaceInt& getFaceInt(unsigned n) const;
  const Polyhedron& getPolyhedron(unsigned n) const;
  unsigned getFacesExtNo() const;
  unsigned getFacesIntNo() const;
  unsigned getPolyhedraNo() const;

  // Prints all the informations about the mesh
  void printAll(std::ostream& out = std::cout) const;

  // Prints only the first five elements of each entity of the mesh.
  void printHead(std::ostream& out = std::cout) const;

  virtual ~Mesh() = default;

  // Proxy used to modify the mesh.
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

  // Function that computes the internal faces of the mesh and completes
  // the information about the external ones.
  void computeFaces();

  // Function that computes the bounding box and the diameter of each polyhedron.
  void computePolyInfo();

  void print(unsigned lineNo, std::ostream& out = std::cout) const;
};

}

#endif // _MESH_HPP_
