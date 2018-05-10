#ifndef _FE_FACE_INT_HPP_
#define _FE_FACE_INT_HPP_

#include "FaceInt.hpp"
#include "FeFace.hpp"
#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace PolyDG
{

class FeFaceInt : public FeFace
{
public:
  FeFaceInt(const FaceInt& face, unsigned order, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule2D& triaRule);

  // Default copy-constructor.
  FeFaceInt(const FeFaceInt&) = default;

  // Default move-constructor.
  FeFaceInt(FeFaceInt&&) = default;

  // Function that returns the value of the basis function f computed at the
  // quadrature point p on the side s.
  inline Real getPhi(SideType s, SizeType p, SizeType f) const;

  // Function that returns the gradient of the basis function f computed at the
  // quadrature point p on the side s.
  inline const Eigen::Vector3d& getPhiDer(SideType s, SizeType p, SizeType f) const;

  // Function that returns the id-number of the elements to which the face belongs
  // from outside.
  inline unsigned getElemOut() const;

  // Function that prints the values of the basis functions over the face.
  void printBasis(std::ostream& out) const override;

  // Function that prints the values of the gradient of the basis functions over
  // the face.
  void printBasisDer(std::ostream& out) const override;

  // Default virtual destructor.
  virtual ~FeFaceInt() = default;

private:
  // Auxiliary function that computes the values of the basis functions and
  // their gradient and fills phi_ and phDer_.
  void compute_basis() override;

  // Auxiliary function that, given the s of the face, the quadrature point p
  // and basis function f, returns the index in which the corresponding value
  // is stored in phi_ and phiDer_.
  inline SizeType sub2ind(SideType s, SizeType p, SizeType f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceInt::getPhi(SideType s, SizeType p, SizeType f) const
{
  return phi_[sub2ind(s, p, f)];
}

inline const Eigen::Vector3d& FeFaceInt::getPhiDer(SideType s, SizeType p, SizeType f) const
{
  return phiDer_[sub2ind(s, p, f)];
}

inline unsigned FeFaceInt::getElemOut() const
{
  return static_cast<const FaceInt&>(face_).getTetOut().getPoly().getId();
}

inline SizeType FeFaceInt::sub2ind(SideType s, SizeType p, SizeType f) const
{
  return (s == In) + 2 * (f + p * dof_);
}

} // namespace PolyDG

#endif // _FE_FACE_INT_HPP_
