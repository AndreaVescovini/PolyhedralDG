#include "FaceExt.hpp"

namespace PolyDG
{

FaceExt::FaceExt(Vertex& v1, Vertex& v2,  Vertex& v3,
                 Real area, Eigen::Vector3d normal, BCtype bcLabel,
                 Tetrahedron& tet1, unsigned faceNoTet1)
  :  FaceExt(v1, v2, v3, area, normal, bcLabel, &tet1, faceNoTet1) {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
                 Real area, Eigen::Vector3d normal, BCtype bcLabel,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, area, normal, tet1, faceNoTet1), bcLabel_{bcLabel} {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCtype bcLabel,
                 Tetrahedron& tet1 , unsigned faceNoTet1)
  :  FaceExt(v1, v2, v3, bcLabel, &tet1, faceNoTet1) {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCtype bcLabel,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, tet1, faceNoTet1), bcLabel_{bcLabel} {}

void FaceExt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId() << ", L: " << bcLabel_
      << ", A:" << areaDoubled_ << ", N:" << normal_.transpose();
}

} // namespace PolyDG
