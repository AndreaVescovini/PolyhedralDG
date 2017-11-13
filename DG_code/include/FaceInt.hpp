#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

namespace geom
{

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

  const Tetrahedron& getTet2() const;
  Tetrahedron& getTet2();
  unsigned getFaceNoTet2() const;
  void setTet2(Tetrahedron& tet2);
  void setTet2(Tetrahedron* tet2);
  void setFaceNoTet2(unsigned faceNoTet2);

private:
  // I store also the second tetrahedron that owns the face and the local number
  // of the face in tetrahderon.
  Tetrahedron* tet2_;
  unsigned faceNoTet2_;

  void print(std::ostream& out) const override;
};

}

#endif // _FACE_INT_HPP_
