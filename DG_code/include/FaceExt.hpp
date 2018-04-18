#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include "PolyDG.hpp"
#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

#include <iostream>

namespace PolyDG
{

class FaceExt : public FaceAbs
{
public:
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
          Real area, Eigen::Vector3d normal, BCtype bcLabel,
          Tetrahedron& tet1, unsigned faceNoTet1);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
          Real area, Eigen::Vector3d normal, BCtype bcLabel,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCtype bcLabel,
          Tetrahedron& tet1 , unsigned faceNoTet1);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, BCtype bcLabel,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  inline BCtype getBClabel() const;
  inline void setBClabel(BCtype bcLabel);

  virtual ~FaceExt() = default;

private:
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
