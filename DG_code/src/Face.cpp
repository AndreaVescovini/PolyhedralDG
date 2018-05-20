/*!
    @file   Face.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class Face
*/

#include "Face.hpp"

#include <Eigen/Geometry>

#include <algorithm>

namespace PolyDG
{

Face::Face(Vertex& v1, Vertex& v2, Vertex& v3)
  : vertices_{v1, v2, v3}, tetIn_{nullptr}, faceNoTetIn_{4}
{
  // I sort vertices comparing the id.
  std::sort(vertices_.begin(), vertices_.end(), compareId);
}

Face::Face(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn, unsigned faceNoTetIn)
  : Face(v1, v2, v3)
{
  tetIn_ = &tetIn;
  faceNoTetIn_ = faceNoTetIn;
}

} // namespace PolyDG
