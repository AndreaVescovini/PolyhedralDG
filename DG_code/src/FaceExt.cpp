#include "FaceExt.hpp"

namespace geom {

FaceExt::FaceExt(const Vertex& v1, const Vertex& v2, const Vertex& v3, unsigned BClabel)
  : Face(v1, v2, v3), BClabel_{BClabel} {}

unsigned FaceExt::getBClabel() const
{
  return BClabel_;
}

void FaceExt::print(std::ostream& out) const
{
  out << id_ << " " << "V: " << vertices_[0].get().getId() << " "
                             << vertices_[1].get().getId() << " "
                             << vertices_[2].get().getId() << ", L: " << BClabel_;
}

}
