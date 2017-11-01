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
  MeshReaderPoly() = default;
  explicit MeshReaderPoly(std::array<std::string, 3> sections);

  void read(Mesh& mesh, const std::string& fileName) const;

  void setSections(std::array<std::string, 3> sections);
  std::array<std::string, 3> getSections() const;

private:
  std::array<std::string, 3> sections_;
};

}

#endif // _MESH_READER_POLY_HPP_
