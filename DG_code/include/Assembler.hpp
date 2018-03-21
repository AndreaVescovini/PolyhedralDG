#ifndef _ASSEMBLER_HPP_
#define _ASSEMBLER_HPP_

#include "FeSpace.hpp"
#include "geom.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "Expr.hpp"

namespace dgfem
{

template <typename T>
class Assembler
{
public:
  Assembler(const Expr<T>& expr, const FeSpace& Vh);

  void assembleVol();

  void printMatrix(std::ostream& out = std::cout) const;

  virtual ~Assembler() = default;

private:
  const Expr<T>& expr_;;
  const FeSpace& Vh_;
  Eigen::SparseMatrix<double> A_;

};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename T>
Assembler<T>::Assembler(const Expr<T>& expr, const FeSpace& Vh)
  : expr_{expr}, Vh_{Vh},
    A_{Vh.getFeElementsNo() * Vh.getDofNo(), Vh.getFeElementsNo() * Vh.getDofNo()} {}

template <typename T>
void Assembler<T>::assembleVol()
{
  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeElementsNo());

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    geom::real sum;

    for(unsigned i = 0; i < Vh_.getDofNo(); i++)
      for(unsigned j = i; j < Vh_.getDofNo(); j++)
      {
        sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += expr_(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + elemNo * Vh_.getDofNo(), j + elemNo * Vh_.getDofNo(), sum);
      }

    elemNo++;
  }

  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numericale zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

template <typename T>
void Assembler<T>::printMatrix(std::ostream& out) const
{
  out << A_.selfadjointView<Eigen::Upper>() << std::endl;
}

}

#endif // _ASSEMBLER_HPP_
