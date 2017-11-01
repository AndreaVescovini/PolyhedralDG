#include "MeshProxy.hpp"

namespace dgfem {

MeshProxy::MeshProxy(Mesh& mesh)
  : mesh_{mesh} {}

std::vector<Point>& MeshProxy::getVerticesRef() const
{
  return mesh_.vertices_;
}

std::vector<Tetrahedron>& MeshProxy::getTetrahedraRef() const
{
  return mesh_.tetrahedra_;
}

}
