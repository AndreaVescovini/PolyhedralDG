/*!
    @file   FaceInt.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class FaceInt
*/

#include "FaceInt.hpp"

namespace PolyDG
{

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3)
  : FaceAbs(v1, v2, v3), tetOut_{nullptr} {}

FaceInt::FaceInt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
                 unsigned faceNoTetIn, Tetrahedron& tetOut)
  : FaceAbs(v1, v2, v3, tetIn, faceNoTetIn), tetOut_{&tetOut} {}

void FaceInt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId()
      << ", T:" << tetIn_->getId() << " " << tetOut_->getId()
      << ", A:" << areaDoubled_ / 2 << ", N:" << normal_.transpose();
}

} // namespace PolyDG
