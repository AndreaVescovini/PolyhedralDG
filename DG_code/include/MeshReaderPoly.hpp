#ifndef _MESH_READER_POLY_HPP_
#define _MESH_READER_POLY_HPP_

#include "MeshReader.hpp"
#include "Mesh.hpp"

#include <string>
#include <array>
#include <fstream>

namespace PolyDG
{

class MeshReaderPoly : public MeshReader
{
public:
  explicit MeshReaderPoly(const std::array<std::string, 4>& sections =
    {"Vertices", "Tetrahedra", "Triangles", "Polyhedra"});

  // Method that reads the file fileName containing the mesh.
  void read(Mesh& mesh, const std::string& fileName) const;

  inline void setSections(const std::array<std::string, 4>& sections);
  inline const std::array<std::string, 4>& getSections() const;

  virtual ~MeshReaderPoly() = default;

private:
  // Sections names in the mesh file.
  std::array<std::string, 4> sections_;

  bool goToSection(std::ifstream& meshFile, unsigned secNo) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline void MeshReaderPoly::setSections(const std::array<std::string, 4>& sections)
{
  sections_ = sections;
}

inline const std::array<std::string, 4>& MeshReaderPoly::getSections() const
{
  return sections_;
}

} // namespace PolyDG

#endif // _MESH_READER_POLY_HPP_
