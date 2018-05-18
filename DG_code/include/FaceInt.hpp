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

    This class defines an internal face of a polyhedral mesh.
    A internal face is defined as one of the co-planar triangles belonging to the
    triangulation of an interface between two polyhedral elements. An interface
    is the intersection of the two-dimensional facets of neighbouring elements.
    This class inherits from FaceAbs and extends it adding a second Tetrahedron
    that shares it. Note that the outward normal give by getNormal() is that one
    pointing from the Tetrahedron "In" to the Tetrhedron "Out".
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

      @param v1 Vertex.
      @param v2 Vertex.
      @param v3 Vertex.
      @param tetIn Tetrahedron In to which the face belongs.
      @param faceNoTetIn Number of the face in the tetrahedron tetIn.
      @apram tetOut Tetrhedron Out to which the face belongs.
  */
  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, Tetrahedron& tetOut);

  //! Copy constructor
  FaceInt(const FaceInt&) = default;

  //! Move constructor
  FaceInt(FaceInt&&) = default;

  /*!
      @brief  Check if a Tetrahedron "Out" is set
      @return @b true if a Tetrahedron "Out" is set, @b false if it is not.
  */
  inline bool isTetOutSet() const;

  /*!
      @brief   Get the Tetrahedron "Out" to which the face belongs
      @warning This function can be called only if isTetOutSet() returns @b true.
  */
  inline const Tetrahedron& getTetOut() const;

  /*!
      @brief   Get the Tetrahedron "Out" to which the face belongs
      @warning This function can be called only if isTetOutSet() returns @b true.
  */
  inline Tetrahedron& getTetOut();

  /*!
      @brief Set a Tetrahedron "Out"
      @param tetIn The Tetrahedron "Out" to which this face belongs.
  */
  inline void setTetOut(Tetrahedron& tetOut);

  //! Destructor
  virtual ~FaceInt() = default;

private:
  //! Pointer to the tetrahedron Out owning the face
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
