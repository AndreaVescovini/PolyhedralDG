#include "FaceExt.hpp"

namespace geom {

FaceExt::FaceExt(const std::array<labelType, 3> vertices, unsigned BClabel)
  : Face(vertices), BClabel_{BClabel} {}

unsigned FaceExt::getBClabel() const
{
  return BClabel_;
}

void FaceExt::print(std::ostream& out) const
{
  out << "V: " << vertices_[0] << " " << vertices_[1] << " " << vertices_[2]
      << ", L: " << BClabel_;
}

}
