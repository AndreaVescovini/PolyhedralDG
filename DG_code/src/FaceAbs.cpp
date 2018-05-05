#include "FaceAbs.hpp"

#include <Eigen/Geometry>

#include <algorithm>

namespace PolyDG
{

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3)
  : Face(v1, v2, v3), id_{counter_}
{
  computeNormalandArea();
  counter_++;
}

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn, unsigned faceNoTetIn)
  : Face(v1, v2, v3, tetIn, faceNoTetIn), id_{counter_}
{
  computeNormalandArea();
  counter_++;
}

void FaceAbs::computeNormalandArea()
{
  Eigen::Vector3d tmp1(vertices_[0].get().getCoords() - vertices_[1].get().getCoords());
  Eigen::Vector3d tmp2(vertices_[2].get().getCoords() - vertices_[1].get().getCoords());

  normal_ = tmp1.cross(tmp2);
  areaDoubled_ = normal_.norm(); // this is the area doubled
  normal_ /= areaDoubled_;
  checkNormalSign();
}

void FaceAbs::checkNormalSign()
{
  // I check the correctness of the sign of the normal vector performing a dot
  // product with a vector pointing to the forth vertex of the tetrahedron, i.e.
  // that one that does not belong to the face. If the dot product is positive
  // mean that the two vectors are in the same direction so I have to revert the
  // sign of the normal.

  if(tetIn_ != nullptr)
  {
    const Vertex& fourthVertex = tetIn_->getVertex(3 - faceNoTetIn_);
    if(normal_.dot(fourthVertex.getCoords() - vertices_[0].get().getCoords()) > 0 )
      normal_ *= -1.0;
  }
}

std::ostream& operator<<(std::ostream& out, const FaceAbs& FaceAbs)
{
  FaceAbs.print(out);
  return out;
}

unsigned FaceAbs::counter_ = 0;

} // namespace PolyDG
