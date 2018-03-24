#ifndef _FACE_EXT_HPP_
#define _FACE_EXT_HPP_

#include <iostream>
#include "FaceAbs.hpp"
#include "Tetrahedron.hpp"

namespace geom {

class FaceExt : public FaceAbs
{
public:
  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
          real area, Eigen::Vector3d normal, unsigned BClabel,
          Tetrahedron& tet1, unsigned faceNoTet1);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3,
          real area, Eigen::Vector3d normal, unsigned BClabel,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, unsigned BClabel,
          Tetrahedron& tet1 , unsigned faceNoTet1);

  FaceExt(Vertex& v1, Vertex& v2, Vertex& v3, unsigned BClabel,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  unsigned getBClabel() const;
  void setBClabel(unsigned BClabel);

  virtual ~FaceExt() = default;

private:
  // I store an unsigned correspoing to the boundary condition assigned to this
  // external face.
  unsigned BClabel_; // mettere l'enum

  void print(std::ostream& out) const override;
};

}

#endif // _FACE_EXT_HPP_
