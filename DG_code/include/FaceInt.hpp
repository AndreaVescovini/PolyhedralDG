#ifndef _FACE_INT_HPP_
#define _FACE_INT_HPP_

#include "Face.hpp"

namespace geom
{

class FaceInt : public Face
{
public:
  // FaceInt() = default;
  FaceInt(const Vertex& v1, const Vertex& v2, const Vertex& v3);

  const Tetrahedron& getTet2() const;
  unsigned getFaceNoTet2() const;
  void setTet2(const Tetrahedron& tet2);
  void setTet2(const Tetrahedron* tet2);
  void setFaceNoTet2(unsigned faceNoTet2);

private:
  Tetrahedron const* tet2_;
  unsigned faceNoTet2_;

  void print(std::ostream& out) const override;

};

}

#endif // _FACE_INT_HPP_
