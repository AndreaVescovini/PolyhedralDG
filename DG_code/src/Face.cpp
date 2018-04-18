#include "Face.hpp"

#include <Eigen/Geometry>

#include <algorithm>

namespace PolyDG {

Face::Face(Vertex& v1, Vertex& v2, Vertex& v3,
           Tetrahedron& tet1, unsigned faceNoTet1)
  : Face(v1, v2, v3, &tet1, faceNoTet1) {}

Face::Face(Vertex& v1, Vertex& v2, Vertex& v3,
           Tetrahedron* tet1, unsigned faceNoTet1)
  : vertices_{{v1, v2, v3}}, tet1_{tet1}, faceNoTet1_{faceNoTet1}
{
  // I sort vertices comparing the id.
  std::sort(vertices_.begin(), vertices_.end(), compId);
}

} // namespace PolyDG
