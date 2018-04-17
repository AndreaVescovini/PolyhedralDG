#ifndef _FACE_ABS_HPP_
#define _FACE_ABS_HPP_

#include <iostream>
#include "Face.hpp"

namespace geom {

class FaceAbs : public Face
{
public:
  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3, real areaDoubled, Eigen::Vector3d normal,
          Tetrahedron& tet1, unsigned faceNoTet1);

  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3, real areaDoubled, Eigen::Vector3d normal,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
          Tetrahedron& tet1, unsigned faceNoTet1);

  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3,
          Tetrahedron* tet1 = nullptr, unsigned faceNoTet1 = 0); // il defoult di faceNo non è bello, zero andrebbe bene se le cose fossero numerate da 1 a N

  // Mancano i copy constructor e i move constructor.

  inline real getAreaDoubled() const;
  inline const Eigen::Vector3d& getNormal() const;
  inline void setAreaDoubled(real areaDoubled);
  inline void setNormal(const Eigen::Vector3d& normal);
  
  // Function that checks weather the normal vector has the right sign or has
  // to be reverted.
  // **tet1_ should be different from nullptr in order to use  this function.**
  void checkNormalSign();

  inline unsigned getId() const;
  static void resetCounter(unsigned counter = 0);

  virtual ~FaceAbs() = default;

  friend std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

protected:
  const unsigned id_;
  real areaDoubled_;
  Eigen::Vector3d normal_;
  static unsigned counter_;

  virtual void print(std::ostream& out) const = 0;
};

std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline real FaceAbs::getAreaDoubled() const
{
  return areaDoubled_;
}

inline const Eigen::Vector3d& FaceAbs::getNormal() const
{
  return normal_;
}

inline void FaceAbs::setAreaDoubled(real areaDoubled)
{
  areaDoubled_ = areaDoubled;
}

inline void FaceAbs::setNormal(const Eigen::Vector3d& normal)
{
  normal_ = normal;
}

inline unsigned FaceAbs::getId() const
{
  return id_;
}

}

#endif // _FACE_ABS_HPP_
