#ifndef _FACE_HPP_
#define _FACE_HPP_

#include "PolyDG.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <functional>
#include <memory>

namespace PolyDG
{

class Face
{
public:
  Face(Vertex& v1, Vertex& v2, Vertex& v3);

  Face(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn , unsigned faceNoTetIn);

  // Default copy constructor.
  Face(const Face&) = default;

  // Deleted copy assignment operator.
  Face& operator=(const Face&) = delete;

  // Default move constructor.
  Face(Face&&) = default;

  // Deleted move assignment operator.
  Face& operator=(Face&&) = delete;

  inline const Vertex& getVertex(SizeType i) const;
  inline Vertex& getVertex(SizeType i);

  inline const Tetrahedron& getTetIn() const;
  inline Tetrahedron& getTetIn();
  inline void setTetIn(Tetrahedron& tetIn);

  inline unsigned getFaceNoTetIn() const;
  inline void setFaceNoTetIn(unsigned faceNoTetIn);

  virtual ~Face() = default;

  inline friend bool operator==(const Face& lhs, const Face& rhs);

protected:
  // Vertices are stored sorted on the id.
  std::array<std::reference_wrapper<Vertex>, 3> vertices_;

  // Pointer to the tetrahedron owning the face.
  Tetrahedron* tetIn_;

  // Local number of the face in the tetrhedron pointed by tetIn_.
  // Locally the i-th face is the face made all the vartices of the tetrahedron
  // but the i-th, i = 0,...,3.
  unsigned faceNoTetIn_;
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
    return (*lhs == *rhs);
  }
};

template<>
struct hash<unique_ptr<PolyDG::Face>>
{
  size_t operator()(const unique_ptr<PolyDG::Face>& f) const
  {
    return (hash<PolyDG::Vertex>()(f->getVertex(0)) +
            hash<PolyDG::Vertex>()(f->getVertex(1)) * 37 +
            hash<PolyDG::Vertex>()(f->getVertex(2)) * 37 * 37 + 23);
  }
};

} // namespace std

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

namespace PolyDG
{

inline const Vertex& Face::getVertex(SizeType i) const
{
  return vertices_[i];
}

inline Vertex& Face::getVertex(SizeType i)
{
  return vertices_[i];
}

inline const Tetrahedron& Face::getTetIn() const
{
  return *tetIn_;
}

inline Tetrahedron& Face::getTetIn()
{
  return *tetIn_;
}

inline void Face::setTetIn(Tetrahedron& tetIn)
{
  tetIn_ = &tetIn;
}

inline unsigned Face::getFaceNoTetIn() const
{
  return faceNoTetIn_;
}

inline void Face::setFaceNoTetIn(unsigned faceNoTetIn)
{
  faceNoTetIn_ = faceNoTetIn;
}

inline bool operator==(const Face& lhs, const Face& rhs)
{
  return (lhs.vertices_[0] == rhs.vertices_[0] &&
          lhs.vertices_[1] == rhs.vertices_[1] &&
          lhs.vertices_[2] == rhs.vertices_[2]);
}

} // namespace PolyDG

#endif // _FACE_HPP_
