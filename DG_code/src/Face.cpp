#include "Face.hpp"

namespace geom {

// Face::Face()
//   : id_{counter_} // devo inizializare anche gli altri memnbri?
//   {
//     counter_++;
//   }

Face::Face(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  : id_{counter_}, vertices_{{v1, v2, v3}}
  {
    counter_++;
    // fare area e normal vector
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
