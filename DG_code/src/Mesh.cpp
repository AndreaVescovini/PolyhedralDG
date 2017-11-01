#include "Mesh.hpp"

#include <algorithm>
#include "MeshReader.hpp"


namespace dgfem {

// void Mesh::load(const std::string& fileName)
// {
//   // With MeshReader I read vertices and tetrahedra data
//   MeshReader mr("../meshes/", "Vertices", "Tetrahedra", "Polyhedra");
//   mr.read(*this, "cube_str48h.mesh");
//
//
// }

Mesh::Mesh(const std::string& fileName, MeshReader& reader)
{
  reader.read(*this, fileName);
}


void Mesh::print(unsigned lineNo, std::ostream& out) const
{
  out << "-------- MESH --------" << std::endl;

  out << "VERTICES: " << vertices_.size() << std::endl;
  for(unsigned i = 0; i < std::min<size_t>(lineNo, vertices_.size()); i++)
    out << vertices_[i] << std::endl;

  out << "\nTETRAHEDRA: " << tetrahedra_.size() << ", POLYHEDRA: " << polyhedra_.size()
      << std::endl;
  for(unsigned i = 0; i < std::min<size_t>(lineNo, tetrahedra_.size()); i++)
    out << tetrahedra_[i] << std::endl;

  out << "----------------------" << std::endl;
}

void Mesh::printAll(std::ostream& out) const
{
  unsigned lineNo = std::max(vertices_.size(), tetrahedra_.size());
  this->print(lineNo, out);
}

void Mesh::printHead(std::ostream& out) const
{
  this->print(5, out);
}

// void Mesh::setVertices(const std::vector<Point>& vertices)
// {
//   vertices_ = vertices;
// }
//
// void Mesh::setTetrahedra(const std::vector<Tetrahedron>& tetrahedra)
// {
//   tetrahedra_ = tetrahedra;
// }
//
// void Mesh::setPolyhedra(const std::vector<Polyhedron>& polyhedra)
// {
//   polyhedra_ = polyhedra;
// }
//
// void Mesh::setVerticesNo(unsigned verticesNo)
// {
//   verticesNo_ = verticesNo;
// }
//
// void Mesh::setTetrahedraNo(unsigned tetrahedraNo)
// {
//   tetrahedraNo_ = tetrahedraNo;
// }
//
// void Mesh::setPolyhedraNo(unsigned polyhedraNo)
// {
//   polyhedraNo_ = polyhedraNo;
// }

}
