#ifndef _FACE_HPP_
#define _FACE_HPP_

#include <array>
#include <iostream>
#include <functional>
#include <Eigen/Dense>
#include "Vertex.hpp"
#include "Tetrahedron.hpp"

namespace geom
{

class Face
{
public:
  // Face();
  Face(const Vertex& v1, const Vertex& v2, const Vertex& v3); // aggiungere defoult arguments

  const Tetrahedron& getTet1() const;
  unsigned getFaceNoTet1() const;
  void setTet1(const Tetrahedron& tet1);
  void setTet1(const Tetrahedron* tet1);
  void setFaceNoTet1(unsigned faceNoTet1);

  unsigned getId() const;

  real getArea() const;
  const Eigen::Vector3d& getNormal() const;

  static void resetCounter(unsigned counter = 0);

  friend std::ostream& operator<<(std::ostream& out, const Face& face);

protected:
  const unsigned id_;
  std::array<std::reference_wrapper<const Vertex>, 3> vertices_;
  Tetrahedron const* tet1_;
  unsigned faceNoTet1_;
  real area_;
  Eigen::Vector3d normal_;

  static unsigned counter_;


  virtual void print(std::ostream& out) const = 0;

};

std::ostream& operator<<(std::ostream& out, const Face& face);

}

#endif // _FACE_HPP_
