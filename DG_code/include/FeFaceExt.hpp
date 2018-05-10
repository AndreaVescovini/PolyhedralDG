#ifndef _FE_FACE_EXT_HPP_
#define _FE_FACE_EXT_HPP_

#include "FaceExt.hpp"
#include "FeFace.hpp"
#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace PolyDG
{

class FeFaceExt : public FeFace
{
public:
  FeFaceExt(const FaceExt& face, unsigned order, unsigned dof,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule2D& triaRule);

  // Default copy-constructor.
  FeFaceExt(const FeFaceExt&) = default;

  // Default move-constructor.
  FeFaceExt(FeFaceExt&&) = default;

  // Function that returns the value of the basis function f computed at the
  // quadrature point p.
  inline Real getPhi(SizeType p, SizeType f) const;

  // Function that returns the gradient of the basis function f computed at the
  // quadrature point p.
  inline const Eigen::Vector3d& getPhiDer(SizeType p, SizeType f) const;

  // Function that returns the type of boundary condition of the face.
  inline BCType getBClabel() const;

  // Function that prints the values of the basis functions over the face.
  void printBasis(std::ostream& out) const override;

  // Function that prints the values of the gradient of the basis functions over
  // the face.
  void printBasisDer(std::ostream& out) const override;

  // Default virtual destructor.
  virtual ~FeFaceExt() = default;

private:
  // Auxiliary function that computes the values of the basis functions and
  // their gradient and fills phi_ and phDer_.
  void compute_basis() override;

  // Auxiliary function that, given the quadrature point p and basis function f,
  // returns the index in which the corresponding value is stored in phi_ and phiDer_.
  inline SizeType sub2ind(SizeType p, SizeType f) const;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real FeFaceExt::getPhi(SizeType p, SizeType f) const
{
  return phi_[sub2ind(p, f)];
}

inline const Eigen::Vector3d& FeFaceExt::getPhiDer(SizeType p, SizeType f) const
{
  return phiDer_[sub2ind(p, f)];
}

inline BCType FeFaceExt::getBClabel() const
{
  return static_cast<const FaceExt&>(face_).getBClabel();
}

inline SizeType FeFaceExt::sub2ind(SizeType p, SizeType f) const
{
  return f + dof_ * p;
}

} // namespace PolyDG

#endif // _FE_FACE_EXT_HPP_
