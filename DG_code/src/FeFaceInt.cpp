#include "FeFaceInt.hpp"
#include "Legendre.hpp"
#include <cmath>

namespace dgfem
{

FeFaceInt::FeFaceInt(const TheFace& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule)
  : FeFace(order, dofNo, basisComposition, triaRule), face_{face}
{
  compute_basis();
}

void FeFaceInt::compute_basis()
{
  Eigen::Vector3d hb = face_.getTet1().getPoly().getBoundingBox().sizes() / 2;
  Eigen::Vector3d mb = face_.getTet1().getPoly().getBoundingBox().center();

  Eigen::Vector3d hb2 = face_.getTet2().getPoly().getBoundingBox().sizes() / 2;
  Eigen::Vector3d mb2 = face_.getTet2().getPoly().getBoundingBox().center();

  unsigned quadPointsNo = triaRule_.getPointsNo();

  phi_.reserve(2 * quadPointsNo * dofNo_);
  phiDer_.reserve(2 * quadPointsNo * dofNo_);

  // loop over quadrature points
  for(unsigned p = 0; p < quadPointsNo; p++)
  {
    // I map the point from the reference tetrahedron to the physical one,
    // then I rescale it in order to compute the scaled legendre polynomial.
    Eigen::Vector3d physicPt = ((face_.getTet1().getMap() *
                                 QuadRuleManager::getFaceMap(face_.getFaceNoTet1()) *
                                 triaRule_.getPoint(p).homogeneous()) - mb).array() / hb.array();
    Eigen::Vector3d physicPt2 = ((face_.getTet1().getMap() *
                                  QuadRuleManager::getFaceMap(face_.getFaceNoTet1()) *
                                  triaRule_.getPoint(p).homogeneous()) - mb2).array() / hb2.array();

    // loop over basis functions
    for(unsigned f = 0; f < dofNo_; f++)
    {
      std::array<std::array<geom::real, 2>, 3> polval;
      std::array<std::array<geom::real, 2>, 3> polval2;

      for(unsigned i = 0; i < 3; i++)
      {
        // I compute the basis for both the sides of the face
        polval[i] = basis::legendre(basisComposition_[f][i], physicPt[i]);
        polval2[i] = basis::legendre(basisComposition_[f][i], physicPt2[i]);

        // I rescale the results in order to have the scaled Legendre polynomials
        polval[i][0] /= std::sqrt(hb[i]);
        polval[i][1] /= (std::sqrt(hb[i]) * hb[i]);
        polval2[i][0] /= std::sqrt(hb2[i]);
        polval2[i][1] /= (std::sqrt(hb2[i]) * hb2[i]);
      }

      phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);
      phi_.emplace_back(polval2[0][0] * polval2[1][0] * polval2[2][0]);

      phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                           polval[0][0] * polval[1][1] * polval[2][0],
                           polval[0][0] * polval[1][0] * polval[2][1]);
      phiDer_.emplace_back(polval2[0][1] * polval2[1][0] * polval2[2][0],
                           polval2[0][0] * polval2[1][1] * polval2[2][0],
                           polval2[0][0] * polval2[1][0] * polval2[2][1]);
    }
  }
}

geom::real FeFaceInt::getPhi(unsigned side, unsigned p, unsigned f) const
{
  return phi_[sub2ind(side, p, f)];
}

const Eigen::Vector3d& FeFaceInt::getPhiDer(unsigned side, unsigned p, unsigned f) const
{
  return phiDer_[sub2ind(side, p, f)];
}

unsigned FeFaceInt::sub2ind(unsigned side, unsigned p, unsigned f) const
{
  return side + 2 * (f + p * dofNo_);
}

void FeFaceInt::printBasis(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << '\n';
  out << face_.getVertex(0).getCoords().transpose() << " - " << face_.getVertex(1).getCoords().transpose() << " - " << face_.getVertex(2).getCoords().transpose() << '\n';

  for(unsigned side = 0; side < 2; side++)
  {
    for(unsigned f = 0; f < dofNo_; f++)
    {
      for(unsigned p = 0; p < triaRule_.getPointsNo(); p++)
        out << getPhi(side, p, f) << ' ';
      out << '\n';
    }
    out << '\n';
  }
  out << '\n';
}

void FeFaceInt::printBasisDer(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << '\n';
  out << face_.getVertex(0).getCoords().transpose() << " - " << face_.getVertex(1).getCoords().transpose() << " - " << face_.getVertex(2).getCoords().transpose() << '\n';

  for(unsigned side = 0; side < 2; side++)
  {
    for(unsigned f = 0; f < dofNo_; f++)
    {
      for(unsigned p = 0; p < triaRule_.getPointsNo(); p++)
        out << getPhiDer(side, p, f).transpose() << '\n';
      out << '\n';
    }
    out << '\n';
  }
  out << '\n';
}

}