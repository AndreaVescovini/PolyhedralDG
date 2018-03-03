#ifndef _FACE_HPP_
#define _FACE_HPP_

#include <array>
#include <functional>
#include <Eigen/Dense>
#include <memory>
#include "geom.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"

namespace geom {

class Face
{
public:
  Face(Vertex& v1, Vertex& v2, Vertex& v3,
       Tetrahedron& tet1 , unsigned faceNoTet1);

  Face(Vertex& v1, Vertex& v2, Vertex& v3,
       Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  // Ci mancano i copy constructor e i move constructor

  const Tetrahedron& getTet1() const;
  Tetrahedron& getTet1();
  unsigned getFaceNoTet1() const;
  void setTet1(Tetrahedron& tet1);
  void setTet1(Tetrahedron* tet1);
  void setFaceNoTet1(unsigned faceNoTet1);

  const Vertex& getVertex(unsigned i) const;
  Vertex& getVertex(unsigned i);

  virtual ~Face() = default;

protected:
  // Vertices are stored sorted on the id.
  std::array<std::reference_wrapper<Vertex>, 3> vertices_;
  // Pointer to the tetrahedron owning the face.
  Tetrahedron* tet1_;
  // Local number of the face in the tetrhedron pointed by tet1_.
  // Locally the i-th face is the face made all the vartices of the tetrahedron
  // but the i-th, i = 0,...,3.
  unsigned faceNoTet1_;
};

}

namespace std {

// I overload the struct equal_to and hash in order to use unordered set of
// unique pointers to Face; I compare faces by the vertices, i.e. by the id
// of vertices.

template<>
struct equal_to<unique_ptr<geom::Face>>
{
  bool operator()(const unique_ptr<geom::Face>& lhs, const unique_ptr<geom::Face>& rhs) const
  {
    bool res = equal_to<geom::Vertex>()(lhs->getVertex(0), rhs->getVertex(0)) &&
               equal_to<geom::Vertex>()(lhs->getVertex(1), rhs->getVertex(1)) &&
               equal_to<geom::Vertex>()(lhs->getVertex(2), rhs->getVertex(2));
    return res;
  }
};

template<>
struct hash<unique_ptr<geom::Face>>
{
  std::size_t operator()(const unique_ptr<geom::Face>& f) const
  {
    std::size_t res = hash<geom::Vertex>()(f->getVertex(0)) +
                      hash<geom::Vertex>()(f->getVertex(1)) * 37 +
                      hash<geom::Vertex>()(f->getVertex(2)) * 37 * 37 + 23;
    return res;
  }
};

}

#endif // _FACE_HPP_
