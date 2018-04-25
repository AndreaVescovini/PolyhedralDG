#include "FaceAbs.hpp"

#include <Eigen/Geometry>

#include <algorithm>

namespace PolyDG
{

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 Real areaDoubled, Eigen::Vector3d normal,
                 Tetrahedron& tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, areaDoubled, normal, &tet1, faceNoTet1) {}

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 Real areaDoubled, Eigen::Vector3d normal,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  Face(v1, v2, v3, tet1, faceNoTet1),
     id_{counter_}, areaDoubled_{areaDoubled}, normal_{normal}
{
  counter_++;
}

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 Tetrahedron& tet1 , unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, &tet1, faceNoTet1) {}

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  Face(v1, v2, v3, tet1, faceNoTet1), id_{counter_}
{
  // I compute the normal vector and the area exploiting the cross product
  // between vectors of edges.
  Eigen::Vector3d tmp1(v1.getCoords() - v2.getCoords());
  Eigen::Vector3d tmp2(v3.getCoords() - v2.getCoords());

  normal_ = tmp1.cross(tmp2);
  areaDoubled_ = normal_.norm(); // this is the area doubled
  normal_ /= areaDoubled_;
  if(tet1_ != nullptr)
    checkNormalSign();

  counter_++;
}

void FaceAbs::checkNormalSign()
{
  // I check the correctness of the sign of the normal vector performing a dot
  // product with a vector pointing to the forth vertex of the tetrahedron, i.e.
  // that one that does not belong to the face. If the dot product is positive
  // mean that the two vectors are in the same direction so I have to revert the
  // sign of the normal.
  const Vertex& fourthVertex = tet1_->getVertex(3 - faceNoTet1_);

  if(normal_.dot(fourthVertex.getCoords() - vertices_[0].get().getCoords()) > 0 )
    normal_ *= -1.0;

  // std::cout << "Checked normal sign" << std::endl;
}

std::ostream& operator<<(std::ostream& out, const FaceAbs& FaceAbs)
{
  FaceAbs.print(out);
  return out;
}

unsigned FaceAbs::counter_ = 0;

} // namespace PolyDG
