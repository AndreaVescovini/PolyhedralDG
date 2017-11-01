#ifndef _MESH_PROXY_HPP_
#define _MESH_PROXY_HPP_

#include <vector>
#include "Mesh.hpp"
#include "Point.hpp"
#include "Tetrahedron.hpp"

namespace dgfem {

class Mesh;

class MeshProxy
{
public:
  explicit MeshProxy(Mesh& mesh);

  std::vector<geom::Point>& getVerticesRef() const;
  std::vector<geom::Tetrahedron>& getTetrahedraRef() const;

private:
  Mesh& mesh_;
};

}

#endif // _MESH_PROXY_HPP_
