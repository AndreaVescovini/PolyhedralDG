#include "FeFaceInt.hpp"
#include "Legendre.hpp"

#include <algorithm>
#include <cmath>

namespace PolyDG
{

FeFaceInt::FeFaceInt(const FaceInt& face, unsigned order, unsigned dof,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule2D& triaRule)
  : FeFace(face, dof, basisComposition, triaRule)
{
  penaltyParam_ = order * order / std::min(face.getTetIn().getPoly().getDiameter(),
                                           face.getTetOut().getPoly().getDiameter());
  compute_basis();
}

void FeFaceInt::compute_basis()
{
  // hb and hb2 contain the half of the dimensions of the bounding box of the two
  // polyhedra sharing the face. They are needed for the computation of the scaled Legendre polynomials.
  const Eigen::Vector3d hbIn  = face_.getTetIn().getPoly().getBoundingBox().sizes() / 2;
  const Eigen::Vector3d hbOut = static_cast<const FaceInt&>(face_).getTetOut().getPoly().getBoundingBox().sizes() / 2;

  // mb and mb2 contain the center of the bounding box of the two polyhedra sharing
  // the face. They are needed for the computation of the scaled Legendre polynomials.
  const Eigen::Vector3d mbIn  = face_.getTetIn().getPoly().getBoundingBox().center();
  const Eigen::Vector3d mbOut = static_cast<const FaceInt&>(face_).getTetOut().getPoly().getBoundingBox().center();

  const SizeType quadPointsNo = triaRule_.getPointsNo();

  phi_.reserve(2 * quadPointsNo * dof_);
  phiDer_.reserve(2 * quadPointsNo * dof_);

  // loop over quadrature points
  for(SizeType p = 0; p < quadPointsNo; p++)
  {
    // I map the quadrature point from the refrence triangle to the face of the
    // reference tetrahedron and then to the physical one, finally I rescale it
    // in order to compute the scaled legendre polynomial.
    const Eigen::Vector3d physicPt = this->getQuadPoint(p);
    const Eigen::Vector3d physicPtIn  = (physicPt - mbIn).array() / hbIn.array();
    const Eigen::Vector3d physicPtOut = (physicPt - mbOut).array() / hbOut.array();

    // Loop over basis functions.
    for(unsigned f = 0; f < dof_; f++)
    {
      std::array<std::array<Real, 2>, 3> polvalIn;
      std::array<std::array<Real, 2>, 3> polvalOut;

      // Loop over the three coordinates.
      for(unsigned i = 0; i < 3; i++)
      {
        // I compute the basis for both the sides of the face.
        polvalIn[i][0] = (legendre(basisComposition_[f][i], physicPtIn(i)) / std::sqrt(hbIn(i)));
        polvalIn[i][1] = (legendreDer(basisComposition_[f][i], physicPtIn(i)) / std::sqrt(hbIn(i)) / hbIn(i));
        polvalOut[i][0] = (legendre(basisComposition_[f][i], physicPtOut(i)) / std::sqrt(hbOut(i)));
        polvalOut[i][1] = (legendreDer(basisComposition_[f][i], physicPtOut(i)) / std::sqrt(hbOut(i)) / hbOut(i));
      }

      phi_.emplace_back(polvalIn[0][0] * polvalIn[1][0] * polvalIn[2][0]);
      phi_.emplace_back(polvalOut[0][0] * polvalOut[1][0] * polvalOut[2][0]);

      phiDer_.emplace_back(polvalIn[0][1] * polvalIn[1][0] * polvalIn[2][0],
                           polvalIn[0][0] * polvalIn[1][1] * polvalIn[2][0],
                           polvalIn[0][0] * polvalIn[1][0] * polvalIn[2][1]);
      phiDer_.emplace_back(polvalOut[0][1] * polvalOut[1][0] * polvalOut[2][0],
                           polvalOut[0][0] * polvalOut[1][1] * polvalOut[2][0],
                           polvalOut[0][0] * polvalOut[1][0] * polvalOut[2][1]);
    }
  }
}

void FeFaceInt::printBasis(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << ": [ ";
  out << face_.getVertex(0).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(1).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(2).getCoords().transpose() << "]\n";

  std::array<SideType, 2> sides = { Out, In };

  // Loop over sides.
  for(SizeType i = 0; i < 2; i++)
  {
    out << "Side " << (i == 0 ? "Out" : "In") << '\n';
    out << "Basis Functions = 1, 2, ..., dof per element" << '\n';
    // Loop over quadrature points.
    for(SizeType p = 0; p < triaRule_.getPointsNo(); p++)
    {
      out << "Quad. Point " << p + 1 << ": ";

      // Loop over basis functions.
      for(unsigned f = 0; f < dof_; f++)
        out << getPhi(sides[i], p, f) << ' ';

      out << '\n';
    }
    out << '\n';
  }
  out << '\n';
}

void FeFaceInt::printBasisDer(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << ": [ ";
  out << face_.getVertex(0).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(1).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(2).getCoords().transpose() << "]\n";

  std::array<SideType, 2> sides = { Out, In };

  // Loop over sides.
  for(SizeType i = 0; i < 2; i++)
  {
    out << "Side " << (i == 0 ? "Out" : "In") << '\n';
    out << "Basis Functions = 1, 2, ..., dof per element" << '\n';
    // Loop over quadrature points.
    for(SizeType p = 0; p < triaRule_.getPointsNo(); p++)
    {
      out << "Quad. Point " << p + 1 << ": ";

      // Loop over basis functions.
      for(unsigned f = 0; f < dof_; f++)
        out << "[ " << getPhiDer(sides[i], p, f).transpose() << " ] ";

      out << '\n';
    }
    out << '\n';
  }
  out << '\n';
}

} // namespace PolyDG
