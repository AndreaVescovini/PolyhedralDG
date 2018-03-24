#ifndef _MESH_PROXY_HPP_
#define _MESH_PROXY_HPP_

#include <vector>
#include "Mesh.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"
#include "FaceExt.hpp"
#include "Polyhedron.hpp"

namespace geom {

// using geom::Vertex;
// using geom::Tetrahedron;
// using geom::FaceExt;
// using geom::Polyhedron;

class Mesh;

class MeshProxy
{
public:
  explicit MeshProxy(Mesh& mesh);

  std::vector<Vertex>& getVerticesRef() const;
  std::vector<Tetrahedron>& getTetrahedraRef() const;
  std::vector<FaceExt>& getFacesExtRef() const;
  std::vector<Polyhedron>& getPolyhedraRef() const;

  virtual ~MeshProxy() = default;

private:
  Mesh& mesh_;
};

}

#endif // _MESH_PROXY_HPP_
