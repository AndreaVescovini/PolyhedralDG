#ifndef _ASSEMBLER_HPP_
#define _ASSEMBLER_HPP_

#include "FeSpace.hpp"
#include "geom.hpp"
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "ExprWrapper.hpp"
#include <iterator>

namespace dgfem
{

class Assembler
{
public:
  explicit Assembler(const FeSpace& Vh);

  template <typename T>
  void assembleVol(const ExprWrapper<T>& expr);

  template <typename T>
  void assembleFacesInt(const ExprWrapper<T>& expr);

  void printMatrix(std::ostream& out = std::cout) const;
  void printMatrixSym(std::ostream& out = std::cout) const;

  virtual ~Assembler() = default;

private:
  const FeSpace& Vh_;
  Eigen::SparseMatrix<double> A_;

};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename T>
void Assembler::assembleVol(const ExprWrapper<T>& expr)
{
// I exploit the conversion to derived
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeElementsNo());

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i <= j; i++)
      {
        geom::real sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + elemNo * Vh_.getDofNo(), j + elemNo * Vh_.getDofNo(), sum);
      }

    elemNo++;
  }

  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

template <typename T>
void Assembler::assembleFacesInt(const ExprWrapper<T>& expr)
{
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeFacesIntNo() * 4);

  for(auto it = Vh_.feFacesIntCbegin(); it != Vh_.feFacesIntCbegin()+1; it++)
    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < Vh_.getDofNo(); i++)//= j; i++)
        for(unsigned side1 = 0; side1 < 2; side1++)
          for(unsigned side2 = 0; side2 < 2; side2++)
          {
            geom::real sum = 0.0;

            for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
              sum += exprDerived(*it, i, j, side1, side2, q) * it->getWeight(q) * it->getAreaDoubled();

              tripletList.emplace_back(i + it->getElem(side1) * Vh_.getDofNo(),
                                       j + it->getElem(side2) * Vh_.getDofNo(),
                                       sum);
              std::cout << sum << '\n';
          }

  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

}

#endif // _ASSEMBLER_HPP_
