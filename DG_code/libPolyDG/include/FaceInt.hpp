/*!
    @file   FaceInt.hpp
    @author Andrea Vescovini
    @brief  Class for internal faces of a polyhedral mesh
*/

#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

#include <iostream>

namespace PolyDG
{

/*!
    @brief Class for internal faces of a polyhedral mesh

    This class defines an internal face of a polyhedral mesh.@n
    A internal face is defined as one of the co-planar triangles belonging to the
    triangulation of an interface between two polyhedral elements. An interface
    is the intersection of the two-dimensional facets of neighbouring elements.@n
    This class inherits from FaceAbs and extends it adding a second Tetrahedron
    that shares it. Note that the outward normal give by getNormal() is that one
    pointing from the Tetrahedron @a "In" to the Tetrhedron @a "Out".
*/

class FaceInt : public FaceAbs
{
public:
  /*!
      @brief Constructor that takes three vertices

      This constructor calls the constructor of FaceAbs.
  */
  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3);

  /*!
      @brief Constructor that takes three vertices two Tetrahedron

      This constructor calls the constructor of FaceAbs and sets the Out
      Tetrhedron..

      @param v1          Vertex.
      @param v2          Vertex.
      @param v3          Vertex.
      @param tetIn       Tetrahedron In to which the face belongs.
      @param faceNoTetIn Number of the face in the tetrahedron tetIn.
      @param tetOut      Tetrhedron Out to which the face belongs.
  */
  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, Tetrahedron& tetOut);

  //! Copy constructor
  FaceInt(const FaceInt&) = default;

  //! Move constructor
  FaceInt(FaceInt&&) = default;

  /*!
      @brief  Check if a Tetrahedron @a "Out" is set
      @return @c true if a Tetrahedron @a "Out" is set, @c false if it is not.
  */
  inline bool isTetOutSet() const;

  /*!
      @brief   Get the Tetrahedron @a "Out" to which the face belongs
      @warning This function can be called only if isTetOutSet() returns @c true.
  */
  inline const Tetrahedron& getTetOut() const;

  /*!
      @brief   Get the Tetrahedron @a "Out" to which the face belongs
      @warning This function can be called only if isTetOutSet() returns @c true.
  */
  inline Tetrahedron& getTetOut();

  /*!
      @brief Set a Tetrahedron @a "Out"
      @param tetOut The Tetrahedron @a "Out" to which this face belongs.
  */
  inline void setTetOut(Tetrahedron& tetOut);

  //! Destructor
  virtual ~FaceInt() = default;

private:
  //! Pointer to the Tetrahedron "Out" owning the face
  Tetrahedron* tetOut_;

  //! Print information about the face
  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline bool FaceInt::isTetOutSet() const
{
  return tetOut_ == nullptr ? false : true;
}

inline const Tetrahedron& FaceInt::getTetOut() const
{
  return *tetOut_;
}

inline Tetrahedron& FaceInt::getTetOut()
{
  return *tetOut_;
}

inline void FaceInt::setTetOut(Tetrahedron& tetOut)
{
  tetOut_ = &tetOut;
}

} // namespace PolyDG

#endif // _FACE_INT_HPP_
