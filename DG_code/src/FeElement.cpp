#include "FeElement.hpp"
#include "Legendre.hpp"
#include <cmath>

namespace dgfem
{

FeElement::FeElement(const theElem& elem, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule<Eigen::Vector3d>& tetraRule)
  : elem_{elem}, order_{order}, dofNo_{dofNo}, tetraRule_{tetraRule}
{
  Eigen::Vector3d hb = elem_.getBoundingBox().sizes() / 2;
  Eigen::Vector3d mb = elem_.getBoundingBox().center();

  unsigned quadPointsNo = tetraRule_.getPointsNo();
  unsigned tetraNo = elem_.getTetrahedraNo();

  phi_.reserve(tetraNo * quadPointsNo * dofNo_);
  phiDer_.reserve(tetraNo * quadPointsNo * dofNo_);

  // loop over tetrahedra
  for(unsigned t = 0; t < tetraNo; t++)
  {
    // loop over quadrature points
    for(unsigned p = 0; p < quadPointsNo; p++)
    {
      // I map the point from the reference tetrahedron to the physical one,
      // then I rescale it in order to compute the scaled legendre polynomial.
      Eigen::Vector3d physicPt = ((elem_.getTetra(t).getMap() * tetraRule_.getPoint(p)) - mb).array() / hb.array();

      // loop over quadrature points
      for(unsigned f = 0; f < dofNo_; f++)
      {
        std::array<std::array<geom::real, 2>, 3> polval;
        for(unsigned i = 0; i < 3; i++)
        {
          polval[i] = basis::legendre(basisComposition[f][i], physicPt[i]);
          polval[i][0] /= sqrt(hb[i]);
          polval[i][1] /= std::sqrt(hb[i]);
          polval[i][1] /= hb[i];
        }

        phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);

        phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                             polval[0][0] * polval[1][1] * polval[2][0],
                             polval[0][0] * polval[1][0] * polval[2][1]);
      }
    }
  }
}

unsigned FeElement::getOrder() const
{
  return order_;
}

unsigned FeElement::getDofNo() const
{
  return dofNo_;
}

geom::real FeElement::getPhi(unsigned t, unsigned p, unsigned f) const
{
  return phi_[sub2ind(t, p, f)];
}

const Eigen::Vector3d& FeElement::getPhiDer(unsigned t, unsigned p, unsigned f) const
{
  return phiDer_[sub2ind(t, p, f)];
}

unsigned FeElement::sub2ind(unsigned t, unsigned p, unsigned f) const
{
  return f + dofNo_ * (p + t * tetraRule_.getPointsNo());
}

}
