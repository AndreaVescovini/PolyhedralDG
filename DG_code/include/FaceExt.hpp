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
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCtype bcLabel);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn,
          unsigned faceNoTetIn, BCtype bcLabel);

  // Default copy constructor.
  FaceExt(const FaceExt&) = default;

  // Default move constructor.
  FaceExt(FaceExt&&) = default;

  inline BCtype getBClabel() const;
  inline void setBClabel(BCtype bcLabel);

  virtual ~FaceExt() = default;

private:
  // Label for a specific boundary condition.
  BCtype bcLabel_;

  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline BCtype FaceExt::getBClabel() const
{
  return bcLabel_;
}

inline void FaceExt::setBClabel(BCtype bcLabel)
{
  bcLabel_ = bcLabel;
}

} // namespace PolyDG

#endif // _FACE_EXT_HPP_
