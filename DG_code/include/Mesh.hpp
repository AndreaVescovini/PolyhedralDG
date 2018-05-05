#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "FaceExt.hpp"
#include "FaceInt.hpp"
#include "MeshProxy.hpp"
#include "MeshReader.hpp"
#include "PolyDG.hpp"
#include "Polyhedron.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace PolyDG
{

class MeshReader;

class Mesh
{
public:
  // Ho tolto i ldefault constructor perchè ho bisogno di chiamare computeFaces
  // e computePolyInfo dopo la lettura del meshFile.
  Mesh(const std::string& fileName, MeshReader& reader);

  Mesh(const Mesh&) = default;
  Mesh& operator=(const Mesh&) = default;
  Mesh(Mesh&&) = default;
  Mesh& operator=(Mesh&&) = default;

  inline const Vertex& getVertex(sizeType n) const;
  inline const Tetrahedron& getTetrahedron(sizeType n) const;
  inline const FaceExt& getFaceExt(sizeType n) const;
  inline const FaceInt& getFaceInt(sizeType n) const;
  inline const Polyhedron& getPolyhedron(sizeType n) const;

  inline sizeType getVerticesNo() const;
  inline sizeType getTetrahedraNo() const;
  inline sizeType getFacesExtNo() const;
  inline sizeType getFacesIntNo() const;
  inline sizeType getPolyhedraNo() const;

  inline Real getMaxDiameter() const;
  inline Real getMinDiameter() const;

  // Prints all the informations about the mesh
  void printAll(std::ostream& out = std::cout) const;

  // Prints only the first five elements of each entity of the mesh.
  void printHead(std::ostream& out = std::cout) const;

  // Prints a piece of information about the mesh.
  void printInfo(std::ostream& out = std::cout) const;

  // Default virtual destructor.
  virtual ~Mesh() = default;

  // Proxy used to modify the mesh.
  friend class MeshProxy;

private:
  std::vector<Vertex> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  std::vector<FaceExt> facesExt_;
  std::vector<FaceInt> facesInt_;
  std::vector<Polyhedron> polyhedra_;

  Real hmax_;
  Real hmin_;

  // Function that computes the internal faces of the mesh and completes
  // the information about the external ones.
  void computeFaces();

  // Function that computes the bounding box and the diameter of each polyhedron.
  void computePolyInfo();

  void print(sizeType lineNo, std::ostream& out = std::cout) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const Vertex& Mesh::getVertex(sizeType n) const
{
  return vertices_[n];
}

inline const Tetrahedron& Mesh::getTetrahedron(sizeType n) const
{
  return tetrahedra_[n];
}

inline const FaceExt& Mesh::getFaceExt(sizeType n) const
{
  return facesExt_[n];
}

inline const FaceInt& Mesh::getFaceInt(sizeType n) const
{
  return facesInt_[n];
}

inline const Polyhedron& Mesh::getPolyhedron(sizeType n) const
{
  return polyhedra_[n];
}

inline sizeType Mesh::getVerticesNo() const
{
  return vertices_.size();
}

inline sizeType Mesh::getTetrahedraNo() const
{
  return tetrahedra_.size();
}

inline sizeType Mesh::getFacesExtNo() const
{
  return facesExt_.size();
}

inline sizeType Mesh::getFacesIntNo() const
{
  return facesInt_.size();
}

inline sizeType Mesh::getPolyhedraNo() const
{
  return polyhedra_.size();
}

inline Real Mesh::getMaxDiameter() const
{
  return hmax_;
}

inline Real Mesh::getMinDiameter() const
{
  return hmin_;
}

} // namespace PolyDG

#endif // _MESH_HPP_
