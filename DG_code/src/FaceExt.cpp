#include "FaceExt.hpp"

namespace PolyDG
{

FaceExt::FaceExt(Vertex& v1, Vertex& v2,  Vertex& v3, BCType bcLabel)
  :  FaceAbs(v1, v2, v3), bcLabel_{bcLabel} {}

FaceExt::FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
                 unsigned faceNoTetIn, BCType bcLabel)
  :  FaceAbs(v1, v2, v3, tetIn, faceNoTetIn), bcLabel_{bcLabel} {}

void FaceExt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId() << ", L: " << bcLabel_
      << ", A:" << areaDoubled_ / 2 << ", N:" << normal_.transpose();
}

} // namespace PolyDG
