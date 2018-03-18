#include "FeFaceExt.hpp"
#include "Legendre.hpp"
#include <cmath>

namespace dgfem
{

FeFaceExt::FeFaceExt(const TheFace& face, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector2d>& triaRule)
  : FeFace(order, dofNo, basisComposition, triaRule), face_{face}
{
  compute_basis();
}

void FeFaceExt::compute_basis()
{
  Eigen::Vector3d hb = face_.getTet1().getPoly().getBoundingBox().sizes() / 2;
  Eigen::Vector3d mb = face_.getTet1().getPoly().getBoundingBox().center();

  unsigned quadPointsNo = triaRule_.getPointsNo();

  phi_.reserve(quadPointsNo * dofNo_);
  phiDer_.reserve(quadPointsNo * dofNo_);

  // loop over quadrature points
  for(unsigned p = 0; p < quadPointsNo; p++)
  {
    // I map the point from the reference tetrahedron to the physical one,
    // then I rescale it in order to compute the scaled legendre polynomial.
    Eigen::Vector3d physicPt = ((face_.getTet1().getMap() *
                                 QuadRuleManager::getFaceMap(face_.getFaceNoTet1()) *
                                 triaRule_.getPoint(p).homogeneous()) - mb).array() / hb.array();

    // loop over basis functions
    for(unsigned f = 0; f < dofNo_; f++)
    {
      std::array<std::array<geom::real, 2>, 3> polval;
      for(unsigned i = 0; i < 3; i++)
      {
        polval[i] = basis::legendre(basisComposition_[f][i], physicPt[i]);

        // I rescale the results in order to have the scaled Legendre polynomials
        polval[i][0] /= std::sqrt(hb[i]);
        polval[i][1] /= (std::sqrt(hb[i]) * hb[i]);
      }

      phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);

      phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                           polval[0][0] * polval[1][1] * polval[2][0],
                           polval[0][0] * polval[1][0] * polval[2][1]);
    }
  }
}

geom::real FeFaceExt::getPhi(unsigned p, unsigned f) const
{
  return phi_[sub2ind(p, f)];
}

const Eigen::Vector3d& FeFaceExt::getPhiDer(unsigned p, unsigned f) const
{
  return phiDer_[sub2ind(p, f)];
}

unsigned FeFaceExt::sub2ind(unsigned p, unsigned f) const
{
  return f + dofNo_ * p;
}

void FeFaceExt::printBasis(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << '\n';
  out << face_.getVertex(0).getCoords().transpose() << " - " << face_.getVertex(1).getCoords().transpose() << " - " << face_.getVertex(2).getCoords().transpose() << '\n';

  for(unsigned f = 0; f < dofNo_; f++)
  {
    for(unsigned p = 0; p < triaRule_.getPointsNo(); p++)
      out << getPhi(p, f) << ' ';
    out << '\n';
  }
  out << '\n';
}

void FeFaceExt::printBasisDer(std::ostream& out = std::cout) const
{
  out << "Face " << face_.getId() << '\n';
  out << face_.getVertex(0).getCoords().transpose() << " - " << face_.getVertex(1).getCoords().transpose() << " - " << face_.getVertex(2).getCoords().transpose() << '\n';

  for(unsigned f = 0; f < dofNo_; f++)
  {
    for(unsigned p = 0; p < triaRule_.getPointsNo(); p++)
      out << getPhiDer(p, f).transpose() << '\n';
    out << '\n';
  }
  out << '\n';
}

}