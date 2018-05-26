/*!
    @file   MeshReaderPoly.hpp
    @author Andrea Vescovini
    @brief  Class that defines a reader for meshes *.mesh
*/

#ifndef _MESH_READER_POLY_HPP_
#define _MESH_READER_POLY_HPP_

#include "Mesh.hpp"
#include "MeshReader.hpp"

#include <array>
#include <fstream>
#include <string>

namespace PolyDG
{

/*!
    @brief Class that defines a reader for meshes *.mesh

    This class inherits from MeshReader and implements a read function that reads
    tridimensional tetrahedral MEDIT text meshes with the extension *.mesh
    (see https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf for more
    information).@n
    In order to provide a polyhedral mesh you should have installed METIS (see
    http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) and you can
    rely on the bash script provided with PolyDG for customizing the mesh
    generating polyhedra.@n
    If you provide to MeshReaderPoly a simple tetraheldral mesh so that the last
    section "Polyhedra" is not found, it is assumed that every tetrahedron is
    also a polyhedron.
*/

class MeshReaderPoly : public MeshReader
{
public:
  /*!
      @brief Constructor that takes the names of the sections of the mesh file

      The first section must describe vertices, the second section tetrahedra,
      the third external faces and tha last polyhedra. The last section is the
      only one that is not compulsory.

      @param sections std::array with the names of the sections of the file.
                      If unspecified they are "Vertices", "Tetrahedra", "Triangles"
                      and "Polyhedra".
  */
  explicit MeshReaderPoly(const std::array<std::string, 4>& sections =
    {"Vertices", "Tetrahedra", "Triangles", "Polyhedra"});

  //! Copy constructor
  MeshReaderPoly(const MeshReaderPoly&) = default;

  //! Copy-assignment operator
  MeshReaderPoly& operator=(MeshReaderPoly&) = default;

  //! Move constructor
  MeshReaderPoly(MeshReaderPoly&&) = default;

  //! Move-assigment operator
  MeshReaderPoly& operator=(MeshReaderPoly&&) = default;

  /*!
      @brief Read the file with the mesh

      This function performs the actual reading from the file fileName and
      through the proxy saves the data in mesh.

      @param mesh     The Mesh you want to fill.
      @param fileName The name of the file that contains the mesh.

      @attention If you provide a wrong file name a @c std::runtime_error
                 exception is thrown.
      @attention If the file format is different from the expected one a
                 MeshFormatError exception is thrown.
  */
  void read(Mesh& mesh, const std::string& fileName) const;

  //! Set the names of the sections of the mesh file
  inline void setSections(const std::array<std::string, 4>& sections);

  //! Get the names of the sections of the mesh file
  inline const std::array<std::string, 4>& getSections() const;

  //! Destructor
  virtual ~MeshReaderPoly() = default;

private:
  //! Sections names in the mesh file
  std::array<std::string, 4> sections_;

  //! Jump to the section secNo while reading the mesh file
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
