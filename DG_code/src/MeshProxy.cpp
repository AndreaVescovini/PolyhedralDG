#include "MeshProxy.hpp"

namespace PolyDG
{

MeshProxy::MeshProxy(Mesh& mesh)
  : mesh_{mesh} {}

std::vector<Vertex>& MeshProxy::getVerticesRef() const
{
  return mesh_.vertices_;
}

std::vector<Tetrahedron>& MeshProxy::getTetrahedraRef() const
{
  return mesh_.tetrahedra_;
}

std::vector<FaceExt>& MeshProxy::getFacesExtRef() const
{
  return mesh_.facesExt_;
}

std::vector<Polyhedron>& MeshProxy::getPolyhedraRef() const
{
  return mesh_.polyhedra_;
}

} // namespace PolyDG
