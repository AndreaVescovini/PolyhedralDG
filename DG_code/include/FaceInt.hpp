#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include <iostream>
#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

namespace geom {

class FaceInt : public FaceAbs
{
public:
  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
          real area, Eigen::Vector3d normal,
          Tetrahedron& tet1, unsigned faceNoTet1,
          Tetrahedron& tet2, unsigned faceNoTet2);

  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
          real area, Eigen::Vector3d normal,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0,
          Tetrahedron* tet2 = nullptr, unsigned faceNoTet2 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
          Tetrahedron& tet1 , unsigned faceNoTet1,
          Tetrahedron& tet2 , unsigned faceNoTet2);

  FaceInt(Vertex& v1, Vertex& v2, Vertex& v3,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0, // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N
          Tetrahedron* tet2 = nullptr, unsigned faceNoTet2 = 0);

  inline const Tetrahedron& getTet2() const;
  inline Tetrahedron& getTet2();
  inline unsigned getFaceNoTet2() const;
  inline void setTet2(Tetrahedron& tet2);
  inline void setTet2(Tetrahedron* tet2);
  inline void setFaceNoTet2(unsigned faceNoTet2);

  virtual ~FaceInt() = default;

private:
  // I store also the second tetrahedron that owns the face and the local number
  // of the face in tetrahderon.
  Tetrahedron* tet2_;
  unsigned faceNoTet2_;

  void print(std::ostream& out) const override;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const Tetrahedron& FaceInt::getTet2() const
{
  return *tet2_;
}

inline Tetrahedron& FaceInt::getTet2()
{
  return *tet2_;
}

inline unsigned FaceInt::getFaceNoTet2() const
{
  return faceNoTet2_;
}

inline void FaceInt::setTet2(Tetrahedron& tet2)
{
  setTet2(&tet2);
}

inline void FaceInt::setTet2(Tetrahedron* tet2)
{
  tet2_ = tet2;
}

inline void FaceInt::setFaceNoTet2(unsigned faceNoTet2)
{
  faceNoTet2_ = faceNoTet2;
}

} // namespace geom

#endif // _FACE_INT_HPP_
