#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include "FaceAbs.hpp"
#include "PolyDG.hpp"
#include "Tetrahedron.hpp"

#include <iostream>

namespace PolyDG
{

class FaceExt : public FaceAbs
{
public:
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCType bcLabel);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, BCType bcLabel);

  // Default copy constructor.
  FaceExt(const FaceExt&) = default;

  // Default move constructor.
  FaceExt(FaceExt&&) = default;

  inline BCType getBClabel() const;
  inline void setBClabel(BCType bcLabel);

  virtual ~FaceExt() = default;

private:
  // Label for a specific boundary condition.
  BCType bcLabel_;

  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline BCType FaceExt::getBClabel() const
{
  return bcLabel_;
}

inline void FaceExt::setBClabel(BCType bcLabel)
{
  bcLabel_ = bcLabel;
}

} // namespace PolyDG

#endif // _FACE_EXT_HPP_
