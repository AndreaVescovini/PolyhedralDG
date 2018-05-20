/*!
    @file   MeshProxy.hpp
    @author Andrea Vescovini
    @brief  Class that defines a proxy for accessing a Mesh
*/

#ifndef _MESH_PROXY_HPP_
#define _MESH_PROXY_HPP_

#include "FaceExt.hpp"
#include "FaceInt.hpp"
#include "Mesh.hpp"
#include "Polyhedron.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <vector>

namespace PolyDG
{

class Mesh;

/*!
    @brief Class that defines a proxy for accessing a Mesh

    This class defines a proxy for accessing the containers of a Mesh.@n
    It should be used in a MeshReader in order to construct a Mesh.

*/
class MeshProxy
{
public:
  /*!
      @brief Cosntructor
      @param mesh The mesh I want to access.
  */
  explicit MeshProxy(Mesh& mesh);

  //! Copy constructor
  MeshProxy(const MeshProxy&) = default;

  //! Deleted copy-assigment operator
  // I do not want that mesh_ is copy assigned.
  MeshProxy& operator=(const MeshProxy&) = delete;

  //! Move constructor
  MeshProxy(MeshProxy&&) = default;

  //! Deleted move-assigment operator
  // I do not want that mesh_ is move assigned.
  MeshProxy& operator=(MeshProxy&&) = delete;

  //! Get a reference to the vector of vertices
  std::vector<Vertex>& getVerticesRef() const;

  //! Get a reference to the vector of tetrahedra
  std::vector<Tetrahedron>& getTetrahedraRef() const;

  //! Get a reference to the vector of external faces
  std::vector<FaceExt>& getFacesExtRef() const;

  //! Get a reference to the vector of internal faces
  std::vector<FaceInt>& getFacesIntRef() const;

  //! Get a reference to the vector of polyhedra
  std::vector<Polyhedron>& getPolyhedraRef() const;

  //! Destructor
  virtual ~MeshProxy() = default;

private:
  //! Reference to the mesh
  Mesh& mesh_;
};

} // namespace PolyDG

#endif // _MESH_PROXY_HPP_
