#ifndef _ASSEMBLER_HPP_
#define _ASSEMBLER_HPP_

#include "FeSpace.hpp"
#include "geom.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "ExprWrapper.hpp"

namespace dgfem
{

class Assembler
{
public:
  explicit Assembler(const FeSpace& Vh);

  template <typename T>
  void assembleVol(const ExprWrapper<T>& expr);

  void printMatrix(std::ostream& out = std::cout) const;

  virtual ~Assembler() = default;

private:
  const FeSpace& Vh_;
  Eigen::SparseMatrix<double> A_;

};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename T>
void Assembler::assembleVol(const ExprWrapper<T>& expr)
{
  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeElementsNo());

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    geom::real sum;

    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i <= j; i++)
      {
        sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += expr(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + elemNo * Vh_.getDofNo(), j + elemNo * Vh_.getDofNo(), sum);
      }

    elemNo++;
  }

  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

}

#endif // _ASSEMBLER_HPP_
