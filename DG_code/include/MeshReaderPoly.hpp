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
  MeshReaderPoly(std::string folder, std::string titleSec1,
                 std::string titleSec2, std::string titleSec3);

  void read(Mesh& mesh, const std::string& fileName) const;

  void setSections(std::string titleSec1, std::string titleSec2, std::string titleSec3);
  std::array<std::string, 3> getSections() const;

private:
  std::array<std::string, 3> sections_;
};

}

#endif // _MESH_READER_POLY_HPP_
