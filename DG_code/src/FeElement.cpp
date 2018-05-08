#include "FeElement.hpp"
#include "Legendre.hpp"

#include <cmath>

namespace PolyDG
{

FeElement::FeElement(const Element& elem, unsigned dof,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRuleManager::Rule3D& tetraRule)
  : elem_{elem}, dof_{dof}, basisComposition_{basisComposition}, tetraRule_{tetraRule}
{
  compute_basis();
}

void FeElement::compute_basis()
{
  // hb contains the half of the dimensions of the bounding box of the polyhedron,
  // mb contains the center of the bounding box. They are needed for the computation
  // of the scaled Legendre polynomials.
  const Eigen::Vector3d hb = elem_.getBoundingBox().sizes() / 2;
  const Eigen::Vector3d mb = elem_.getBoundingBox().center();

  const SizeType quadPointsNo = tetraRule_.getPointsNo();
  const SizeType tetraNo = elem_.getTetrahedraNo();

  phi_.reserve(tetraNo * quadPointsNo * dof_);
  phiDer_.reserve(tetraNo * quadPointsNo * dof_);

  // Loop over tetrahedra.
  for(SizeType t = 0; t < tetraNo; t++)
  {
    // Loop over quadrature points.
    for(SizeType p = 0; p < quadPointsNo; p++)
    {
      // I map the quadrature point from the reference tetrahedron to the physical one,
      // then I rescale it in order to compute the scaled legendre polynomial.
      const Eigen::Vector3d physicPt = ((elem_.getTetra(t).getMap() * tetraRule_.getPoint(p)) - mb).array() / hb.array();

      // loop over basis functions.
      for(unsigned f = 0; f < dof_; f++)
      {
        std::array<std::array<Real, 2>, 3> polval;

        // Loop over coordinates.
        for(unsigned i = 0; i < 3; i++)
        {
          // polval stores the value of the Legendre polynomial and its derivative.
          polval[i][0] = (legendre(basisComposition_[f][i], physicPt(i)) / std::sqrt(hb(i)));
          polval[i][1] = (legendreDer(basisComposition_[f][i], physicPt(i)) / std::sqrt(hb(i)) / hb(i));
        }

        phi_.emplace_back(polval[0][0] * polval[1][0] * polval[2][0]);

        // The gradient is computed deriving monomials one by one.
        phiDer_.emplace_back(polval[0][1] * polval[1][0] * polval[2][0],
                             polval[0][0] * polval[1][1] * polval[2][0],
                             polval[0][0] * polval[1][0] * polval[2][1]);
      }
    }
  }
}

void FeElement::printBasis(std::ostream& out) const
{
  // Loop over tetrahedra.
  for(SizeType t = 0; t < elem_.getTetrahedraNo(); t++)
  {
    out << "Tetrahedron " << elem_.getTetra(t).getId() << '\n';

    // Loop over quadrature points.
    for(SizeType p = 0; p < tetraRule_.getPointsNo(); p++)
    {
      // Loop over basis functions.
      for(SizeType f = 0; f < dof_; f++)
        out << getPhi(t, p, f) << ' ';

      out << '\n';
    }
    out << '\n';
  }
}

void FeElement::printBasisDer(std::ostream& out) const
{
  // Loop over tetrahedra.
  for(SizeType t = 0; t < elem_.getTetrahedraNo(); t++)
  {
    out << "Tetrahedron " << elem_.getTetra(t).getId() << '\n';

    // Loop over quadrature points.
    for(SizeType p = 0; p < tetraRule_.getPointsNo(); p++)
    {

      // Loop over basis functions.
      for(SizeType f = 0; f < dof_; f++)
        out << getPhiDer(t, p, f).transpose() << '\n';

      out << '\n';
    }
    out << '\n';
  }
}

} // namespace PolyDG
