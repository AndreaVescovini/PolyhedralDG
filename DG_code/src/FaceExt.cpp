#include "FaceExt.hpp"

namespace geom {

FaceExt::FaceExt(Vertex& v1, Vertex& v2,  Vertex& v3,
                 real area, Eigen::Vector3d normal, unsigned BClabel,
                 Tetrahedron& tet1, unsigned faceNoTet1)
  :  FaceExt(v1, v2, v3, area, normal, BClabel, &tet1, faceNoTet1) {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
                 real area, Eigen::Vector3d normal, unsigned BClabel,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, area, normal, tet1, faceNoTet1), BClabel_{BClabel} {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, unsigned BClabel,
                 Tetrahedron& tet1 , unsigned faceNoTet1)
  :  FaceExt(v1, v2, v3, BClabel, &tet1, faceNoTet1) {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, unsigned BClabel,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, tet1, faceNoTet1), BClabel_{BClabel} {}

void FaceExt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId() << ", L: " << BClabel_
      << ", A:" << areaDoubled_ << ", N:" << normal_.transpose();
}

} // namespace geom
