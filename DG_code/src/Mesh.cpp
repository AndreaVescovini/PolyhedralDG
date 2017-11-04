#include "Mesh.hpp"

#include <algorithm>

namespace dgfem {

Mesh::Mesh(const std::string& fileName, MeshReader& reader)
{
  reader.read(*this, fileName);
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

void Mesh::print(unsigned lineNo, std::ostream& out) const
{
  out << "-------- MESH --------" << std::endl;

  out << "VERTICES: " << vertices_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, vertices_.size()); i++)
    out << vertices_[i] << std::endl;

  out << "\nTETRAHEDRA: " << tetrahedra_.size() << ", POLYHEDRA: " << polyhedra_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, tetrahedra_.size()); i++)
    out << tetrahedra_[i] << std::endl;

  out << "\nEXTERNAL FACES: " << facesExt_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, facesExt_.size()); i++)
    out << facesExt_[i] << std::endl;

  out << "----------------------" << std::endl;
}

}
