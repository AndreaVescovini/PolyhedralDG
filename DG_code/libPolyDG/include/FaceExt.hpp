/*!
    @file   FaceExt.hpp
    @author Andrea Vescovini
    @brief  Class for external faces of a polyhedral mesh
*/

#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include "FaceAbs.hpp"
#include "PolyDG.hpp"
#include "Tetrahedron.hpp"

#include <iostream>

namespace PolyDG
{

/*!
    @brief Class for external faces of a polyhedral mesh

    This class defines an external face of a polyhedral mesh.@n
    An external face is defined as one of the co-planar triangles belonging to
    the triangulation of the intersection between a facet of a polyhedral element
    and the external boundary of the domain.@n
    This class inherits from FaceAbs and extends it with a label useful to
    define the boundary condition.
*/

class FaceExt : public FaceAbs
{
public:
  /*!
      @brief Constructor that takes three vertices and a label

      This constructor calls the constructor of FaceAbs and sets the label.

      @param v1      Vertex.
      @param v2      Vertex.
      @param v3      Vertex.
      @param bcLabel The label to be set.
  */
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCLabelType bcLabel);

  /*!
      @brief Constructor that takes three vertices, a Tetrahedron and a label

      This constructor calls the constructor of FaceAbs and sets the label.

      @param v1          Vertex.
      @param v2          Vertex.
      @param v3          Vertex.
      @param tetIn       Tetrahedron @a "In" to which the face belongs.
      @param faceNoTetIn Number of the face in the Tetrahedron tetIn.
      @param bcLabel     The label to be set.
  */
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, BCLabelType bcLabel);

  //! Copy constructor
  FaceExt(const FaceExt&) = default;

  //! Move constructor
  FaceExt(FaceExt&&) = default;

  //! Get the label
  inline BCLabelType getBClabel() const;

  //! Set the label
  inline void setBClabel(BCLabelType bcLabel);

  //! Destructor
  virtual ~FaceExt() = default;

private:
  //! Label for a specific boundary condition
  BCLabelType bcLabel_;

  //! Print information about the face
  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline BCLabelType FaceExt::getBClabel() const
{
  return bcLabel_;
}

inline void FaceExt::setBClabel(BCLabelType bcLabel)
{
  bcLabel_ = bcLabel;
}

} // namespace PolyDG

#endif // _FACE_EXT_HPP_
