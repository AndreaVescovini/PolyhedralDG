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

class MeshProxy
{
public:

  // Cosntructor.
  explicit MeshProxy(Mesh& mesh);

  // Default copy constructor.
  MeshProxy(const MeshProxy&) = default;

  // Deleted copy assigment since I do not want that mesh_ is copy assigned.
  MeshProxy& operator=(const MeshProxy&) = delete;

  // Default move constructor.
  MeshProxy(MeshProxy&&) = default;

  // Deleted move assigment since I do not want that mesh_ is move assigned.
  MeshProxy& operator=(MeshProxy&&) = delete;
  
  std::vector<Vertex>& getVerticesRef() const;
  std::vector<Tetrahedron>& getTetrahedraRef() const;
  std::vector<FaceExt>& getFacesExtRef() const;
  std::vector<FaceInt>& getFacesIntRef() const;
  std::vector<Polyhedron>& getPolyhedraRef() const;

  virtual ~MeshProxy() = default;

private:
  Mesh& mesh_;
};

} // namespace PolyDG

#endif // _MESH_PROXY_HPP_
