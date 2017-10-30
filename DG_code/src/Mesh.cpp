#include "Mesh.hpp"

#include <algorithm>
#include "MeshReader.hpp"


namespace geom {

void Mesh::load(const std::string& fileName)
{
  // With MeshReader I read vertices and tetrahedra data
  MeshReader mr("../meshes/", "Vertices", "Tetrahedra", "Polyhedra");
  mr.read(*this, "cube_str48h.mesh");


}

void Mesh::print(unsigned lineNo, std::ostream& out) const
{
  out << "-------- MESH --------" << std::endl;

  out << "VERTICES: " << verticesNo_ << std::endl;
  for(unsigned i = 0; i < std::min(lineNo, verticesNo_); i++)
    out << vertices_[i] << std::endl;

  out << "\nTETRAHEDRA: " << tetrahedraNo_ << ", POLYHEDRA: " << polyhedraNo_
      << std::endl;
  for(unsigned i = 0; i < std::min(lineNo, tetrahedraNo_); i++)
    out << tetrahedra_[i] << std::endl;

  out << "----------------------" << std::endl;
}

void Mesh::printAll(std::ostream& out) const
{
  unsigned lineNo = std::max(verticesNo_, tetrahedraNo_);
  this->print(lineNo, out);
}

void Mesh::printHead(std::ostream& out) const
{
  this->print(5, out);
}

void Mesh::setVertices(const std::vector<Point>& vertices)
{
  vertices_ = vertices;
}

void Mesh::setTetrahedra(const std::vector<Tetrahedron>& tetrahedra)
{
  tetrahedra_ = tetrahedra;
}

void Mesh::setPolyhedra(const std::vector<Polyhedron>& polyhedra)
{
  polyhedra_ = polyhedra;
}

void Mesh::setVerticesNo(unsigned verticesNo)
{
  verticesNo_ = verticesNo;
}

void Mesh::setTetrahedraNo(unsigned tetrahedraNo)
{
  tetrahedraNo_ = tetrahedraNo;
}

void Mesh::setPolyhedraNo(unsigned polyhedraNo)
{
  polyhedraNo_ = polyhedraNo;
}

}
