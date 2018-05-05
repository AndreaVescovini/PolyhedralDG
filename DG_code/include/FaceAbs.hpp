#ifndef _FACE_ABS_HPP_
#define _FACE_ABS_HPP_

#include "Face.hpp"

#include <iostream>

namespace PolyDG
{

class FaceAbs : public Face
{
public:
  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3);

  FaceAbs(Vertex& v1, Vertex& v2, Vertex& v3, Tetrahedron& tetIn, unsigned faceNoTetIn);

  // Default copy constructor.
  FaceAbs(const FaceAbs&) = default;

  // Default move constructor.
  FaceAbs(FaceAbs&&) = default;

  inline Real getAreaDoubled() const;
  inline const Eigen::Vector3d& getNormal() const;

  // Function that computes the unitary normal vector to the face and the doubled
  // area of the face.
  void computeNormalandArea();

  // Function that checks weather the normal vector has the right sign (outtward
  // wrt tetIn_) and revertes it if needed.
  void checkNormalSign();

  inline unsigned getId() const;
  inline static void resetCounter(unsigned counter = 0);

  virtual ~FaceAbs() = default;

  friend std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

protected:
  // Id of the face.
  const unsigned id_;

  // Doubled area of the face.
  Real areaDoubled_;

  // Unitary vector normal to the face in the outward direction wrt tetIn_.
  Eigen::Vector3d normal_;

  static unsigned counter_;

  virtual void print(std::ostream& out) const = 0;
};

std::ostream& operator<<(std::ostream& out, const FaceAbs& faceAbs);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FaceAbs::getAreaDoubled() const
{
  return areaDoubled_;
}

inline const Eigen::Vector3d& FaceAbs::getNormal() const
{
  return normal_;
}

inline unsigned FaceAbs::getId() const
{
  return id_;
}

inline void FaceAbs::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _FACE_ABS_HPP_
