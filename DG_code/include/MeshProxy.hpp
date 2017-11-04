#ifndef _MESH_PROXY_HPP_
#define _MESH_PROXY_HPP_

#include <vector>
#include "Mesh.hpp"
#include "Point.hpp"
#include "Tetrahedron.hpp"
#include "FaceExt.hpp"

namespace dgfem {

using geom::Point;
using geom::Tetrahedron;
using geom::FaceExt;

class Mesh;

class MeshProxy
{
public:
  explicit MeshProxy(Mesh& mesh);

  std::vector<Point>& getVerticesRef() const;
  std::vector<Tetrahedron>& getTetrahedraRef() const;
  std::vector<FaceExt>& getFacesExtRef() const;

private:
  Mesh& mesh_;
};

}

#endif // _MESH_PROXY_HPP_
