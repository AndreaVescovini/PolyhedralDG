#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

#include <iostream>

namespace PolyDG
{

class FaceInt : public FaceAbs
{
public:
  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3);

  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, Tetrahedron& tetOut);

  // Default copy constructor.
  FaceInt(const FaceInt&) = default;

  // Default move constructor.
  FaceInt(FaceInt&&) = default;

  inline const Tetrahedron& getTetOut() const;
  inline Tetrahedron& getTetOut();
  inline void setTetOut(Tetrahedron& tetOut);

  virtual ~FaceInt() = default;

private:
  // Pointer to the second tetrahedron that owns the face.
  Tetrahedron* tetOut_;

  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

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
