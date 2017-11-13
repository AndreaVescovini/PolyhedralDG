#include "Face.hpp"
#include <algorithm>
#include <Eigen/Geometry>

namespace geom {

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

const Tetrahedron& Face::getTet1() const
{
  return *tet1_;
}

Tetrahedron& Face::getTet1()
{
  return *tet1_;
}

unsigned Face::getFaceNoTet1() const
{
  return faceNoTet1_;
}

void Face::setTet1(Tetrahedron& tet1)
{
  tet1_ = &tet1;
}

void Face::setTet1(Tetrahedron* tet1)
{
  tet1_ = tet1;
}

void Face::setFaceNoTet1(unsigned faceNoTet1)
{
  faceNoTet1_ = faceNoTet1;
}

const Vertex& Face::getVertex(unsigned i) const
{
  return vertices_[i];
}

Vertex& Face::getVertex(unsigned i)
{
  return vertices_[i];
}

}
