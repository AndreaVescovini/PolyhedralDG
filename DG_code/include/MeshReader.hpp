/*!
    @file   MeshReader.hpp
    @author Andrea Vescovini
    @brief  Abstract class for a reader of a mesh file
*/

#ifndef _MESH_READER_HPP_
#define _MESH_READER_HPP_

#include "Mesh.hpp"

#include <string>

namespace PolyDG
{

class Mesh;

/*!
    @brief Abstract class for a reader of a mesh file

    This base class defines a reader of a mesh file. Every reader should derive
    from this class and give an implementation for the method read(Mesh& mesh, const std::string& fileName).
    This method has to read the mesh file and, through the proxy MeshProxy, fill
    the mesh creating vertices, tetrhedra, polyhedra and external faces
    ( FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCType bcLabel) ). If leaved
    empty, internal faces will be computed later automatically.
*/

class MeshReader
{
public:
  //! Constructor
  MeshReader() = default;

  //! Copy constructor
  MeshReader(const MeshReader&) = default;

  //! Copy-assignment operator
  MeshReader& operator=(MeshReader&) = default;

  //! Move constructor
  MeshReader(MeshReader&&) = default;

  //! Move-assigment operator
  MeshReader& operator=(MeshReader&&) = default;

  /*!
      @brief Read the file with the mesh

      This function performs the actual reading from the file fileName and
      through the proxy saves the data in mesh.

      @param mesh     The Mesh you want to fill.
      @param fileName The name of the file that contains the mesh.
  */
  virtual void read(Mesh& mesh, const std::string& fileName) const = 0;

  //! Destructor
  virtual ~MeshReader() = default;
};

} // namespace PolyDG

#endif // _MESH_READER_HPP_
