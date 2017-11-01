#include "MeshProxy.hpp"

namespace dgfem {

MeshProxy::MeshProxy(Mesh& mesh)
  : mesh_{mesh} {}

std::vector<geom::Point>& MeshProxy::getVerticesRef() const
{
  return mesh_.vertices_;
}

std::vector<geom::Tetrahedron>& MeshProxy::getTetrahedraRef() const
{
  return mesh_.tetrahedra_;
}

}
