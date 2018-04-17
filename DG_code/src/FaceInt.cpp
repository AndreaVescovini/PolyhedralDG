#include "FaceInt.hpp"

namespace geom {

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
                 real area, Eigen::Vector3d normal,
                 Tetrahedron& tet1, unsigned faceNoTet1,
                 Tetrahedron& tet2, unsigned faceNoTet2)
  :  FaceInt(v1, v2, v3, area, normal, &tet1, faceNoTet1, &tet2, faceNoTet2) {}

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
                 real area, Eigen::Vector3d normal,
                 Tetrahedron* tet1, unsigned faceNoTet1,
                 Tetrahedron* tet2, unsigned faceNoTet2)
  :  FaceAbs(v1, v2, v3, area, normal, tet1, faceNoTet1),
     tet2_{tet2}, faceNoTet2_{faceNoTet2} {}

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
                 Tetrahedron& tet1 , unsigned faceNoTet1,
                 Tetrahedron& tet2 , unsigned faceNoTet2)
  :  FaceInt(v1, v2, v3, &tet1, faceNoTet1, &tet2, faceNoTet2) {}

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
                 Tetrahedron* tet1, unsigned faceNoTet1,
                 Tetrahedron* tet2, unsigned faceNoTet2)
  :  FaceAbs(v1, v2, v3, tet1, faceNoTet1),
     tet2_{tet2}, faceNoTet2_{faceNoTet2} {}

void FaceInt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId()
      << ", T:" << tet1_->getId() << " " << tet2_->getId()
      << ", A:" << areaDoubled_ << ", N:" << normal_.transpose();
}

} // namespace geom
