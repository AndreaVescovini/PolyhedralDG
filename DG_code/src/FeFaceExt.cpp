#include "FeFaceExt.hpp"
#include "Legendre.hpp"

#include <cmath>

namespace PolyDG
{

FeFaceExt::FeFaceExt(const FaceExt& face, unsigned order, unsigned dof,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule2D& triaRule)
  : FeFace(face, dof, basisComposition, triaRule)
{
  penaltyParam_ = order * order / face.getTetIn().getPoly().getDiameter();
  compute_basis();
}

void FeFaceExt::compute_basis()
{
  // hb contains the half of the dimensions of the bounding box of the polyhedron,
  // mb contains the center of the bounding box. They are needed for the computation
  // of the scaled Legendre polynomials.
  const Eigen::Vector3d hb = face_.getTetIn().getPoly().getBoundingBox().sizes() / 2;
  const Eigen::Vector3d mb = face_.getTetIn().getPoly().getBoundingBox().center();

  const SizeType quadPointsNo = triaRule_.getPointsNo();

  phi_.reserve(quadPointsNo * dof_);
  phiDer_.reserve(quadPointsNo * dof_);

  // Loop over quadrature points
  for(SizeType p = 0; p < quadPointsNo; p++)
  {
    // I map the quadrature point from the refrence triangle to the face of the
    // reference tetrahedron and then to the physical one, finally I rescale it
    // in order to compute the scaled legendre polynomial.
    const Eigen::Vector3d physicPt = (this->getQuadPoint(p) - mb).array() / hb.array();

    // Loop over basis functions.
    for(unsigned f = 0; f < dof_; f++)
    {
      std::array<std::array<Real, 2>, 3> polval;
      // Loop over the three coordinates.
      for(unsigned i = 0; i < 3; i++)
      {
        polval[i][0] = (legendre(basisComposition_[f][i], physicPt(i)) / std::sqrt(hb(i)));
        polval[i][1] = (legendreDer(basisComposition_[f][i], physicPt(i)) / std::sqrt(hb(i)) / hb(i));
      }

      phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);

      phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                           polval[0][0] * polval[1][1] * polval[2][0],
                           polval[0][0] * polval[1][0] * polval[2][1]);
    }
  }
}

void FeFaceExt::printBasis(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << ": [ ";
  out << face_.getVertex(0).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(1).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(2).getCoords().transpose() << "]\n";

  out << "Basis Functions = 1, 2, ..., dof per element" << '\n';

  // Loop over the quadrature points.
  for(SizeType p = 0; p < triaRule_.getPointsNo(); p++)
  {
    out << "Quad. Point " << p + 1 << ": ";

    // Loop over the basis functions.
    for(unsigned f = 0; f < dof_; f++)
      out << getPhi(p, f) << ' ';

    out << '\n';
  }
  out << '\n';
}

void FeFaceExt::printBasisDer(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << ": [ ";
  out << face_.getVertex(0).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(1).getCoords().transpose() << " ] [ ";
  out << face_.getVertex(2).getCoords().transpose() << "]\n";

  out << "Basis Functions = 1, 2, ..., dof per element" << '\n';

  // Loop over the quadrature points.
  for(SizeType p = 0; p < triaRule_.getPointsNo(); p++)
  {
    out << "Quad. Point " << p + 1 << ": ";

    // Loop over the basis functions.
    for(unsigned f = 0; f < dof_; f++)
      out << "[ " << getPhiDer(p, f).transpose() << " ] ";

    out << '\n';
  }
  out << '\n';
}

} // namespace PolyDG
