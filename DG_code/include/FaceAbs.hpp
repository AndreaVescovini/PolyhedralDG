/*!
    @file   FaceAbs.hpp
    @author Andrea Vescovini
    @brief  Abstract base class for faces of a polyhedral mesh
*/

#ifndef _FACE_ABS_HPP_
#define _FACE_ABS_HPP_

#include "Face.hpp"

#include <iostream>

namespace PolyDG
{

/*!
    @brief Abstract base class for faces of a polyhedral mesh

    This abstract class is the base class for external faces FaceExt and internal
    faces FaceInt of a polyhedral mesh.@n
    A face is defined as one of the co-planar triangles belonging to the
    triangulation of an interface between two polyhedral elements. An interface
    is the intersection of the two-dimensional facets of neighbouring elements.@n
    This class inherits from Face and extends it adding a id number, the area of
    the face and the unitary normal vector, outward with respect to the
    Tetrhedron @a "In" given by getTetIn().
*/

class FaceAbs : public Face
{
public:
  /*!
      @brief Constructor that takes three vertices

      This constructor calls the constructor of Face and sets the id number,
      the normal vector and the area are initialized. The direction of the
      normal vector is random since a Tetrahedron @a "In" has not been set.
  */
  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3);

  /*!
      @brief Constructor that takes three vertices and a Tetrahedron

      This constructor calls the constructor of Face and sets the id number,
      the normal vector and the area are initialized.

      @param v1          Vertex.
      @param v2          Vertex.
      @param v3          Vertex.
      @param tetIn       Tetrahedron @a "In" to which the face belongs.
      @param faceNoTetIn Number of the face in the Tetrahedron tetIn.
  */
  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn, unsigned faceNoTetIn);

  //! Copy constructor.
  FaceAbs(const FaceAbs&) = default;

  //! Move constructor.
  FaceAbs(FaceAbs&&) = default;

  //! Get the the measure of the area doubled
  inline Real getAreaDoubled() const;

  /*!
      @brief  Get the normal vector
      @return @c Eigen::Vector3d containing the unitary normal vector, outward with
              respect to the Tetrhedron @a "In" given by getTetIn().
  */
  inline const Eigen::Vector3d& getNormal() const;

  /*!
      @brief Compute the normal and the area

      This function computes the unitary normal vector to the face, outward with
      respect to the Tetrhedron @a "In" given by getTetIn(), and the doubled area
      of the face. This function calls checkNormalSign().
  */
  void computeNormalandArea();

  /*!
      @brief Check the sign of the normal vector

      This function checks weather the normal vector has the right sign (outward
      wrt the Tetrahedron @a "In") and revertes it if needed. If a Tetrahedron @a "In"
      has not been set this function does nothing.
  */
  void checkNormalSign();

  //! Get the id number
  inline unsigned getId() const;

  /*!
      @brief Reset the counter for id numbers

      This function resets the counter that assigns id numbers. It should be used
      before reading of a new mesh, to assure that id numbers start from 0.

      @param counter The number to which the counter is reset. If not specified
                     it is 0.
  */
  inline static void resetCounter(unsigned counter = 0);

  //! Destructor
  virtual ~FaceAbs() = default;

  /*!
      @brief Overload for the ostream operator

      It calls a function print that is pure virtual in this class.
  */
  friend std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

protected:
  //! Id number
  const unsigned id_;

  //! Doubled area of the face
  Real areaDoubled_;

  //! Unitary vector normal to the face in the outward direction wrt tetIn_
  Eigen::Vector3d normal_;

  //! Counter used to assign a different id to each FaceAbs when it is created
  static unsigned counter_;

  //! Function that prints information about the FaceAbs, pure virtual
  virtual void print(std::ostream& out) const = 0;
};

std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FaceAbs::getAreaDoubled() const
{
  return areaDoubled_;
}

inline const Eigen::Vector3d& FaceAbs::getNormal() const
{
  return normal_;
}

inline unsigned FaceAbs::getId() const
{
  return id_;
}

inline void FaceAbs::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _FACE_ABS_HPP_
