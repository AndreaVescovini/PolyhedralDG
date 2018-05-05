#ifndef _MESH_READER_POLY_HPP_
#define _MESH_READER_POLY_HPP_

#include "Mesh.hpp"
#include "MeshReader.hpp"

#include <array>
#include <fstream>
#include <string>

namespace PolyDG
{

class MeshReaderPoly : public MeshReader
{
public:

  // Constructor that takes the names of the sections of the mesh file.
  explicit MeshReaderPoly(const std::array<std::string, 4>& sections =
    {"Vertices", "Tetrahedra", "Triangles", "Polyhedra"});

  MeshReaderPoly(const MeshReaderPoly&) = default;
  MeshReaderPoly& operator=(MeshReaderPoly&) = default;
  MeshReaderPoly(MeshReaderPoly&&) = default;
  MeshReaderPoly& operator=(MeshReaderPoly&&) = default;

  // Method that reads the file fileName containing the mesh.
  // Exit with value 1 -> Can't open mesh file
  // Exit with value 2 -> Wrong section in the mesh file
  void read(Mesh& mesh, const std::string& fileName) const;

  // Function that sets the names of the sections of the mesh file.
  inline void setSections(const std::array<std::string, 4>& sections);

  // Function that returns the names of the sections of the mesh file.
  inline const std::array<std::string, 4>& getSections() const;

  virtual ~MeshReaderPoly() = default;

private:
  // Sections names in the mesh file.
  std::array<std::string, 4> sections_;

  // Function that jumps to the section secNo while reading the mesh file.
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
