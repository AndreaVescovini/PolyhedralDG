#ifndef _FACE_HPP_
#define _FACE_HPP_

#include "PolyDG.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"

#include <Eigen/Core>

#include <array>
#include <functional>
#include <memory>

namespace PolyDG
{

class Face
{
public:
  Face(Vertex& v1, Vertex& v2, Vertex& v3,
       Tetrahedron& tet1 , unsigned faceNoTet1);

  Face(Vertex& v1, Vertex& v2, Vertex& v3,
       Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  // Ci mancano i copy constructor e i move constructor

  inline const Tetrahedron& getTet1() const;
  inline Tetrahedron& getTet1();
  inline unsigned getFaceNoTet1() const;
  inline void setTet1(Tetrahedron& tet1);
  inline void setTet1(Tetrahedron* tet1);
  inline void setFaceNoTet1(unsigned faceNoTet1);

  inline const Vertex& getVertex(unsigned i) const;
  inline Vertex& getVertex(unsigned i);

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

} // namespace PolyDG

namespace std {

// I overload the struct equal_to and hash in order to use unordered set of
// unique pointers to Face; I compare faces by the vertices, i.e. by the id
// of vertices.

template<>
struct equal_to<unique_ptr<PolyDG::Face>>
{
  bool operator()(const unique_ptr<PolyDG::Face>& lhs, const unique_ptr<PolyDG::Face>& rhs) const
  {
    bool res = equal_to<PolyDG::Vertex>()(lhs->getVertex(0), rhs->getVertex(0)) &&
               equal_to<PolyDG::Vertex>()(lhs->getVertex(1), rhs->getVertex(1)) &&
               equal_to<PolyDG::Vertex>()(lhs->getVertex(2), rhs->getVertex(2));
    return res;
  }
};

template<>
struct hash<unique_ptr<PolyDG::Face>>
{
  std::size_t operator()(const unique_ptr<PolyDG::Face>& f) const
  {
    std::size_t res = hash<PolyDG::Vertex>()(f->getVertex(0)) +
                      hash<PolyDG::Vertex>()(f->getVertex(1)) * 37 +
                      hash<PolyDG::Vertex>()(f->getVertex(2)) * 37 * 37 + 23;
    return res;
  }
};

} // namespace std

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

namespace PolyDG
{

inline const Tetrahedron& Face::getTet1() const
{
  return *tet1_;
}

inline Tetrahedron& Face::getTet1()
{
  return *tet1_;
}

inline unsigned Face::getFaceNoTet1() const
{
  return faceNoTet1_;
}

inline void Face::setTet1(Tetrahedron& tet1)
{
  tet1_ = &tet1;
}

inline void Face::setTet1(Tetrahedron* tet1)
{
  tet1_ = tet1;
}

inline void Face::setFaceNoTet1(unsigned faceNoTet1)
{
  faceNoTet1_ = faceNoTet1;
}

inline const Vertex& Face::getVertex(unsigned i) const
{
  return vertices_[i];
}

inline Vertex& Face::getVertex(unsigned i)
{
  return vertices_[i];
}

} // namespace PolyDG

#endif // _FACE_HPP_
