#include "FaceInt.hpp"

namespace geom {

FaceInt::FaceInt(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  : Face(v1, v2, v3) {}

const Tetrahedron& FaceInt::getTet2() const
{
  return *tet2_;
}

unsigned FaceInt::getFaceNoTet2() const
{
  return faceNoTet2_;
}

void FaceInt::setTet2(const Tetrahedron& tet2)
{
  tet1_ = &tet2;
}

void FaceInt::setTet2(const Tetrahedron* tet2)
{
  tet2_ = tet2;
}

void FaceInt::setFaceNoTet2(unsigned faceNoTet2)
{
  faceNoTet2_ = faceNoTet2;
}

void FaceInt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId();
}

}
