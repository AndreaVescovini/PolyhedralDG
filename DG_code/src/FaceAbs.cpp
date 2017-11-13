#include "FaceAbs.hpp"
#include <algorithm>
#include <Eigen/Geometry>

namespace geom {

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 real area, Eigen::Vector3d normal,
                 Tetrahedron& tet1, unsigned faceNoTet1)
  :  FaceAbs(v1, v2, v3, area, normal, &tet1, faceNoTet1) {}

FaceAbs::FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
                 real area, Eigen::Vector3d normal,
                 Tetrahedron* tet1, unsigned faceNoTet1)
  :  Face(v1, v2, v3, tet1, faceNoTet1),
     id_{counter_}, area_{area}, normal_{normal}
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
  area_ = normal_.norm(); // this is the area doubled
  normal_ /= area_;
  if(tet1_ != nullptr)
    checkNormalSign();

  counter_++;
}

real FaceAbs::getArea() const
{
  return area_;
}

const Eigen::Vector3d& FaceAbs::getNormal() const
{
  return normal_;
}

void FaceAbs::setArea(real area)
{
  area_ = area;
}

void FaceAbs::setNormal(const Eigen::Vector3d& normal)
{
  normal_ = normal;
}

void FaceAbs::checkNormalSign()
{
  // I check the correctness of the sign of the normal vector performing a dot
  // product with a vector pointing to the forth vertex of the tetrahedron, i.e.
  // that one that does not belong to the face. If the dot product is positive
  // mean that the two vectors are in the same direction so I have to revert the
  // sign of the normal.
  const Vertex& fourthVertex = tet1_->getVertex(faceNoTet1_);

  if(normal_.dot(fourthVertex.getCoords() - vertices_[0].get().getCoords()) > 0 )
    normal_ *= -1;

  std::cout << "Checked normal sign" << std::endl;
}

unsigned FaceAbs::getId() const
{
  return id_;
}

void FaceAbs::resetCounter(unsigned counter)
{
  counter_ = counter;
}

std::ostream& operator<<(std::ostream& out, const FaceAbs& FaceAbs)
{
  FaceAbs.print(out);
  return out;
}

unsigned FaceAbs::counter_ = 0;

}
