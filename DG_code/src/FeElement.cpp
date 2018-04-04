#include "FeElement.hpp"
#include "Legendre.hpp"
#include <cmath>

namespace dgfem
{

FeElement::FeElement(const TheElem& elem, unsigned order, unsigned dofNo,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRuleManager::Rule3D& tetraRule)
  : elem_{elem}, dofNo_{dofNo}, basisComposition_{basisComposition}, tetraRule_{tetraRule}
{
  compute_basis();
}

void FeElement::compute_basis()
{
// hb contains the half of the dimensions of the bounding box of the polyhedron,
// mb contains the center of the bounding box. They are needed for the computation
// of the scaled Legendre polynomials.
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
      // I map the quadrature point from the reference tetrahedron to the physical one,
      // then I rescale it in order to compute the scaled legendre polynomial.
      Eigen::Vector3d physicPt = ((elem_.getTetra(t).getMap() * tetraRule_.getPoint(p)) - mb).array() / hb.array();

      // loop over basis functions
      for(unsigned f = 0; f < dofNo_; f++)
      {
        std::array<std::array<geom::real, 2>, 3> polval;
        for(unsigned i = 0; i < 3; i++)
        {
          // polval stores the value of the Legendre polynomial and its derivative
          polval[i] = basis::legendre(basisComposition_[f][i], physicPt[i]);

          // I rescale the results in order to have the scaled Legendre polynomials
          polval[i][0] /= std::sqrt(hb[i]);
          polval[i][1] /= (std::sqrt(hb[i]) * hb[i]);
        }

        phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);

        // The gradient is computed deriving monomials one by one
        phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                             polval[0][0] * polval[1][1] * polval[2][0],
                             polval[0][0] * polval[1][0] * polval[2][1]);
      }
    }
  }
}

unsigned FeElement::getDofNo() const
{
  return dofNo_;
}

unsigned FeElement::getTetrahedraNo() const
{
  return elem_.getTetrahedraNo();
}

geom::real FeElement::getAbsDetJac(unsigned i) const
{
  return elem_.getTetra(i).getAbsDetJacobian();
}

unsigned FeElement::getQuadPointsNo() const
{
  return tetraRule_.getPointsNo();
}

geom::real FeElement::getWeight(unsigned q) const
{
  return tetraRule_.getWeight(q);
}

const Eigen::Vector3d& FeElement::getQuadPoint(unsigned t, unsigned q) const
{
  return elem_.getTetra(t).getMap() * tetraRule_.getPoint(q);
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

void FeElement::printBasis(std::ostream& out) const
{
  // Loop over tetrahedra
  for(unsigned t = 0; t < elem_.getTetrahedraNo(); t++)
  {
    out << "Tetrahedron " << elem_.getTetra(t).getId() << '\n';

    // Loop over basis functions
    for(unsigned f = 0; f < dofNo_; f++)
    {

      // Loop over quadrature points
      for(unsigned p = 0; p < tetraRule_.getPointsNo(); p++)
        out << getPhi(t, p, f) << ' ';

      out << '\n';
    }
    out << '\n';
  }
}

void FeElement::printBasisDer(std::ostream& out) const
{
  // Loop over tetrahedra
  for(unsigned t = 0; t < elem_.getTetrahedraNo(); t++)
  {
    out << "Tetrahedron " << elem_.getTetra(t).getId() << '\n';

    // Loop over quadrature points
    for(unsigned p = 0; p < tetraRule_.getPointsNo(); p++)
    {

      // Loop over basis functions
      for(unsigned f = 0; f < dofNo_; f++)
        out << getPhiDer(t, p, f).transpose() << '\n';

      out << '\n';
    }
    out << '\n';
  }
}

} // namespace dgfem
