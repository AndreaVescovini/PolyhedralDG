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
  void assembleVol(const ExprWrapper<T>& expr, const bool sym = false);

  template <typename T>
  void assembleFacesInt(const ExprWrapper<T>& expr, const bool sym = false);

  template <typename T>
  void assembleFacesExt(const ExprWrapper<T>& expr, unsigned BClabel = 1, const bool sym = false);

  template <typename T>
  void assembleVolRhs(const ExprWrapper<T>& expr);

  template <typename T>
  void assembleFacesExtRhs(const ExprWrapper<T>& expr, unsigned BClabel = 1);

  void printMatrix(std::ostream& out = std::cout) const;
  void printMatrixSym(std::ostream& out = std::cout) const;
  void printRhs(std::ostream& out = std::cout) const;

  virtual ~Assembler() = default;

private:
  const FeSpace& Vh_;
  Eigen::SparseMatrix<double> A_;
  Eigen::VectorXd b_;

};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename T>
void Assembler::assembleVol(const ExprWrapper<T>& expr, const bool sym)
{
// I exploit the conversion to derived
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getFeElementsNo() * Vh_.getDofNo() * Vh_.getDofNo()); //devo calcolare quanto spazio mi serve

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < (j+1) * sym + (1-sym) * Vh_.getDofNo(); i++)
      {
        geom::real sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + elemNo * Vh_.getDofNo(), j + elemNo * Vh_.getDofNo(), sum);
      }

    elemNo++;
  }

  // probabilmente mi serve fare anche un A_.reserve()
  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

template <typename T>
void Assembler::assembleFacesInt(const ExprWrapper<T>& expr, const bool sym)
{
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeFacesIntNo() * 4); // devo calcolare quanto spazio mi serve

  for(auto it = Vh_.feFacesIntCbegin(); it != Vh_.feFacesIntCend(); it++)
    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < Vh_.getDofNo()*(1-sym)+(j+1)*sym; i++)
        for(int side1 = 0; side1 < 2; side1++)
          for(int side2 = 0; side2 < 2; side2++)
          {
            geom::real sum = 0.0;

            for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
              sum += exprDerived(*it, i, j, side1, side2, q) * it->getWeight(q) * it->getAreaDoubled();

            tripletList.emplace_back(i + it->getElem(side1) * Vh_.getDofNo(),
                                     j + it->getElem(side2) * Vh_.getDofNo(),
                                     sum);
          }

// probabilmente mi serve fare anche un A_.reserve()
  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

template <typename T>
void Assembler::assembleFacesExt(const ExprWrapper<T>& expr, unsigned BClabel, const bool sym)
{
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;
  tripletList.reserve(Vh_.getDofNo() * Vh_.getFeFacesExtNo() * 4);  // devo calcolare quanto spazio mi serve

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == BClabel)
      for(unsigned j = 0; j < Vh_.getDofNo(); j++)
        for(unsigned i = 0; i < Vh_.getDofNo()*(1-sym)+(j+1)*sym; i++)
        {
          geom::real sum = 0.0;

          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, q) * it->getWeight(q) * it->getAreaDoubled();

          tripletList.emplace_back(i + it->getElem() * Vh_.getDofNo(),
                                   j + it->getElem() * Vh_.getDofNo(),
                                   sum);
        }

  A_.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  A_.prune(A_.coeff(0,0));
}

template <typename T>
void Assembler::assembleVolRhs(const ExprWrapper<T>& expr)
{
  const T& exprDerived(expr);

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    for(unsigned i = 0; i < Vh_.getDofNo(); i++)
      for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
          b_(i + elemNo * Vh_.getDofNo()) += exprDerived(*it, i, t, q) *
                                             it->getWeight(q) *
                                             it->getAbsDetJac(t);
    elemNo++;
  }

}

template <typename T>
void Assembler::assembleFacesExtRhs(const ExprWrapper<T>& expr, unsigned BClabel)
{
  const T& exprDerived(expr);

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == BClabel)
      for(unsigned i = 0; i < Vh_.getDofNo(); i++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
          b_(i + it->getElem() * Vh_.getDofNo()) += exprDerived(*it, i, q) *
                                                    it->getWeight(q) *
                                                    it->getAreaDoubled();
}

} // namespace dgfem

#endif // _ASSEMBLER_HPP_
