#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "MeshProxy.hpp"
#include "MeshReader.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"
#include "FaceExt.hpp"
#include "FaceInt.hpp"
#include "Polyhedron.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace PolyDG
{

class MeshReader;

class Mesh
{
public:
  Mesh() = default;

  Mesh(const std::string& fileName, MeshReader& reader);

  inline const FaceExt& getFaceExt(unsigned n) const;
  inline const FaceInt& getFaceInt(unsigned n) const;
  inline const Polyhedron& getPolyhedron(unsigned n) const;

  inline unsigned getFacesExtNo() const;
  inline unsigned getFacesIntNo() const;
  inline unsigned getTetrahedraNo() const;
  inline unsigned getPolyhedraNo() const;

  Real getMaxDiameter() const;

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

  // Function that computes the internal faces of the mesh and completes
  // the information about the external ones.
  void computeFaces();

  // Function that computes the bounding box and the diameter of each polyhedron.
  void computePolyInfo();

  void print(unsigned lineNo, std::ostream& out = std::cout) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const FaceExt& Mesh::getFaceExt(unsigned n) const
{
  return facesExt_[n];
}

inline const FaceInt& Mesh::getFaceInt(unsigned n) const
{
  return facesInt_[n];
}

inline const Polyhedron& Mesh::getPolyhedron(unsigned n) const
{
  return polyhedra_[n];
}

inline unsigned Mesh::getFacesExtNo() const
{
  return facesExt_.size();
}

inline unsigned Mesh::getFacesIntNo() const
{
  return facesInt_.size();
}

inline unsigned Mesh::getTetrahedraNo() const
{
  return tetrahedra_.size();
}

inline unsigned Mesh::getPolyhedraNo() const
{
  return polyhedra_.size();
}

} // namespace PolyDG

#endif // _MESH_HPP_
