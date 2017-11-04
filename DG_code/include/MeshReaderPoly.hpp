#ifndef _MESH_READER_POLY_HPP_
#define _MESH_READER_POLY_HPP_

#include <string>
#include <array>
#include "MeshReader.hpp"
#include "Mesh.hpp"

namespace dgfem {

class MeshReaderPoly : public MeshReader
{
public:
  explicit MeshReaderPoly(std::array<std::string, 4> sections =
    {"Vertices", "Tetrahedra", "Triangles", "Polyhedra"});

  void read(Mesh& mesh, const std::string& fileName) const;

  void setSections(std::array<std::string, 4> sections);
  std::array<std::string, 4> getSections() const;

private:
  std::array<std::string, 4> sections_;

  void goToSection(std::ifstream& meshFile, unsigned secNo) const;
};

}

#endif // _MESH_READER_POLY_HPP_
