#include "Face.hpp"

#include <Eigen/Geometry>

namespace geom {

// Face::Face()
//   : id_{counter_} // devo inizializare anche gli altri memnbri?
//   {
//     counter_++;
//   }

Face::Face(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  : id_{counter_}, vertices_{{v1, v2, v3}}
  {
    Eigen::Vector3d tmp1(v1.getCoords() - v2.getCoords());
    Eigen::Vector3d tmp2(v3.getCoords() - v2.getCoords());
    normal_ = tmp1.cross(tmp2);
    area_ = normal_.norm();
    normal_.normalize();

    counter_++;
  }

const Tetrahedron& Face::getTet1() const
{
  return *tet1_;
}

unsigned Face::getFaceNoTet1() const
{
  return faceNoTet1_;
}

void Face::setTet1(const Tetrahedron& tet1)
{
  tet1_ = &tet1;
}

void Face::setTet1(const Tetrahedron* tet1)
{
  tet1_ = tet1;
}

void Face::setFaceNoTet1(unsigned faceNoTet1)
{
  faceNoTet1_ = faceNoTet1;
}

real Face::getArea() const
{
  return area_;
}

const Eigen::Vector3d& Face::getNormal() const
{
  return normal_;
}

void Face::resetCounter(unsigned counter)
{
  counter_ = counter;
}


std::ostream& operator<<(std::ostream& out, const Face& face)
{
  face.print(out);
  return out;
}

unsigned Face::counter_ = 0;

}
