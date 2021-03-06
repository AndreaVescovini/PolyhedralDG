/*!
    @file   Face.hpp
    @author Andrea Vescovini
    @brief  Class that defines a generic face of a polyhedral mesh
*/

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

/*!
    @brief Class that defines a generic face of a polyhedral mesh

    This base class defines a generic face of a polyhedral mesh.@n
    A face is defined as one of the co-planar triangles belonging to the
    triangulation of an interface between two polyhedral elements. An interface
    is the intersection of the two-dimensional facets of neighbouring elements.@n
    It stores the three vertices of the triangle, a pointer to a Tetrahedron
    to which the face belongs and the number of this face in it. Conventionally
    the i-th face is that one without the (3-i)-th vertex.
    For example the face 1 is made by the vertices 0, 2 and 3.@n
    Vertices are stored sorted by their id number.
*/

class Face
{
public:
  /*!
      @brief Constructor that takes three vertices

      The three vertices are setted sorted by their id number.
  */
  Face(Vertex& v1, Vertex& v2, Vertex& v3);

  /*!
      @brief Constructor that takes three vertices and a Tetrahedron @a "In"

      The three vertices are setted sorted by their id number.

      @param v1          Vertex.
      @param v2          Vertex.
      @param v3          Vertex.
      @param tetIn       Tetrahedron @a "In" to which the face belongs.
      @param faceNoTetIn Number of the face in the Tetrhedron tetIn.
  */
  Face(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn, unsigned faceNoTetIn);

  //! Copy constructor
  Face(const Face&) = default;

  //! Deleted copy-assignment operator
  Face& operator=(const Face&) = delete;

  //! Move constructor
  Face(Face&&) = default;

  //! Deleted move-assignment operator
  Face& operator=(Face&&) = delete;

  /*!
      @brief Get a Vertex

      This functions returns the i-th Vertex.

      @param i The index of the Vertex required, it can be 0, 1 or 2.
  */
  inline const Vertex& getVertex(SizeType i) const;

  /*!
      @brief Get a Vertex

      This functions returns the i-th Vertex.

      @param i The index of the Vertex required, it can be 0, 1 or 2.
  */
  inline Vertex& getVertex(SizeType i);

  /*!
      @brief  Check if a Tetrahedron @a "In" is set
      @return @c true if a Tetrahedron @a "In" is set, @c false if it is not.
  */
  inline bool isTetInSet() const;

  /*!
      @brief   Get the Tetrahedron @a "In" to which the face belongs
      @warning This function can be called only if isTetInSet() returns @c true.
  */
  inline const Tetrahedron& getTetIn() const;

  /*!
      @brief   Get the Tetrahedron @a "In" to which the face belongs
      @warning This function can be called only if isTetInSet() returns @c true.
  */
  inline Tetrahedron& getTetIn();

  /*!
      @brief Set a Tetrahedron @a "In"
      @param tetIn The Tetrahedron @a "In" to which this face belongs.
  */
  inline void setTetIn(Tetrahedron& tetIn);

  /*!
      @brief   Get the number of the face in the Tetrahedron @a "In" to which it belongs

      This functions returns the number of the face in the Tetrahedron @a "In" to which
      it belongs. Conventionally the i-th face is that one without the (3-i)-th
      vertex.

      @return This function should be called only if isTetInSet() returns @c true,
               otherwise it returns a meaningless value like 4.
  */
  inline unsigned getFaceNoTetIn() const;

  /*!
      @brief Set the number of the face in the Tetrahedron @a "In" to which it belongs
      @param faceNoTetIn The number you want to set.
  */
  inline void setFaceNoTetIn(unsigned faceNoTetIn);

  //! Destructor
  virtual ~Face() = default;

  /*!
      @brief  Overload of the operator==
      @return @c true if the three vertices are the same.
  */
  inline friend bool operator==(const Face& lhs, const Face& rhs);

protected:
  //! Vertices of the face, sorted on the id number
  std::array<std::reference_wrapper<Vertex>, 3> vertices_;

  //! Pointer to the Tetrahedron @a "In" owning the face
  Tetrahedron* tetIn_;

  //! Local number of the face in the Tetrhedron @a "In" pointed by tetIn_
  // Locally the i-th face is the face made all the vartices of the tetrahedron
  // but the i-th, i = 0,...,3.
  unsigned faceNoTetIn_;
};

} // namespace PolyDG

namespace std {

// I overload the struct equal_to and hash in order to use unordered set of
// unique pointers to Face; I compare faces by the vertices, i.e. by the id
// of vertices.
//! Specialization of the functor @c std::equal_to for @c std::unique_ptr<PolyDG::Face>
template<>
struct equal_to<unique_ptr<PolyDG::Face>>
{
  //! Call operator
  bool operator()(const unique_ptr<PolyDG::Face>& lhs, const unique_ptr<PolyDG::Face>& rhs) const
  {
    return (*lhs == *rhs);
  }
};

//! Specialization of the functor @c std::hash for @c std::unique_ptr<PolyDG::Face>
template<>
struct hash<unique_ptr<PolyDG::Face>>
{
  //! Call operator
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

inline bool Face::isTetInSet() const
{
  return tetIn_ == nullptr ? false : true;
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
