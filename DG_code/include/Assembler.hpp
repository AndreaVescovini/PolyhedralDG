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
#include <algorithm>

namespace dgfem
{

class Assembler
{
public:
  explicit Assembler(const FeSpace& Vh, bool sym = false);

  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool symExpr = false);

  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool symExpr = false);

  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr, unsigned BClabel = 1, bool symExpr = false);

  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);

  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr, unsigned BClabel = 1);

  void solveLU();
  void solveChol();
  void solveCG(const Eigen::VectorXd& x0, unsigned iterMax,
               geom::real tol = Eigen::NumTraits<geom::real>::epsilon());

  void isSymmetric(bool sym);

  const Eigen::SparseMatrix<geom::real> getMatrix() const;
  const Eigen::VectorXd getRhs() const;
  const Eigen::VectorXd getSolution() const;
  unsigned getDim() const;

  void clearMatrix();
  void clearRhs();

  void printMatrix(std::ostream& out = std::cout) const;
  void printMatrixSym(std::ostream& out = std::cout) const;
  void printRhs(std::ostream& out = std::cout) const;

  virtual ~Assembler() = default;

private:
  const FeSpace& Vh_;
  unsigned dim_;
  Eigen::SparseMatrix<geom::real> A_;
  // std::vector<triplet> tripletList
  Eigen::VectorXd b_;
  Eigen::VectorXd u_;

  bool sym_;


};

////////////////////////////////////////////////////////////////////////////////
//-------------------------------IMPLEMENTATION-------------------------------//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
void Assembler::integrateVol(const ExprWrapper<T>& expr, bool symExpr)
{
  // I exploit the conversion to derived
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;

  // If the variational form is symmetric I store only half elements
  if(symExpr == true)
    tripletList.reserve(dim_ * (Vh_.getDofNo() + 1) / 2);
  else
    tripletList.reserve(dim_ * Vh_.getDofNo());

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = elemNo * Vh_.getDofNo();

    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < (symExpr == true ? j + 1 : Vh_.getDofNo()); i++)
      {
        geom::real sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + indexOffset, j + indexOffset, sum);
      }

    elemNo++;
  }

  Eigen::SparseMatrix<geom::real> curA(dim_, dim_);

  curA.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  curA.prune(A_.coeff(0,0));

  if(sym_ == symExpr)
    A_ += curA;
  else if(sym_ == false && symExpr == true)
    A_ += curA.selfadjointView<Eigen::Upper>();
  else
  {
    std::cerr << "The symmetry has been broken" << std::endl;
    sym_ = false;
    A_ = A_.selfadjointView<Eigen::Upper>();
    A_ += curA;
  }
}

template <typename T>
void Assembler::integrateFacesInt(const ExprWrapper<T>& expr, bool symExpr)
{
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;

  // If the variational form is symmetric I store only half elements
  if(symExpr == true)
    tripletList.reserve(Vh_.getFeFacesIntNo() * Vh_.getDofNo() * (Vh_.getDofNo() + 1) * 2);
  else
    tripletList.reserve(Vh_.getFeFacesIntNo() * Vh_.getDofNo() * Vh_.getDofNo() * 4);

  for(auto it = Vh_.feFacesIntCbegin(); it != Vh_.feFacesIntCend(); it++)
  {
    std::array<unsigned, 2> indexOffset = { it->getElem(0) * Vh_.getDofNo(),
                                            it->getElem(1) * Vh_.getDofNo() };

    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < (symExpr == true ? j + 1 : Vh_.getDofNo()); i++)
        for(int side1 = 0; side1 < 2; side1++)
          for(int side2 = 0; side2 < 2; side2++)
          {
            if(symExpr == true && indexOffset[side1] > indexOffset[side2] && i == j)
              continue;

            geom::real sum = 0.0;

            for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
              sum += exprDerived(*it, i, j, side1, side2, q) * it->getWeight(q) * it->getAreaDoubled();

            unsigned rowIndex = (symExpr == true ? std::min(i + indexOffset[side1], j + indexOffset[side2]) : i + indexOffset[side1] );
            unsigned colIndex = (symExpr == true ? std::max(i + indexOffset[side1], j + indexOffset[side2]) : j + indexOffset[side2] );

            tripletList.emplace_back(rowIndex, colIndex, sum);
          }
  }

  Eigen::SparseMatrix<geom::real> curA(dim_, dim_);

  curA.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  curA.prune(A_.coeff(0,0));

  if(sym_ == symExpr)
    A_ += curA;
  else if(sym_ == false && symExpr == true)
    A_ += curA.selfadjointView<Eigen::Upper>();
  else
  {
    std::cerr << "The symmetry has been broken" << std::endl;
    sym_ = false;
    A_ = A_.selfadjointView<Eigen::Upper>();
    A_ += curA;
  }
}

template <typename T>
void Assembler::integrateFacesExt(const ExprWrapper<T>& expr, unsigned BClabel, bool symExpr)
{
  const T& exprDerived(expr);

  using triplet = Eigen::Triplet<double>;
  std::vector<triplet> tripletList;

  // If the variational form is symmetric I store only half elements,
  // I overestimate considering all the external faces with the same type of
  // boundary conditions
  if(symExpr == true)
    tripletList.reserve(Vh_.getFeFacesExtNo() * Vh_.getDofNo() * (Vh_.getDofNo() + 1) / 2);
  else
    tripletList.reserve(Vh_.getFeFacesExtNo() * Vh_.getDofNo() * Vh_.getDofNo());

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == BClabel)
    {
      unsigned indexOffset = it->getElem() * Vh_.getDofNo();

      for(unsigned j = 0; j < Vh_.getDofNo(); j++)
        for(unsigned i = 0; i < (symExpr == true ? j + 1 : Vh_.getDofNo()); i++)
        {
          geom::real sum = 0.0;

          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, q) * it->getWeight(q) * it->getAreaDoubled();

          tripletList.emplace_back(i + indexOffset, j + indexOffset, sum);
        }
    }

  // I made and overestimate so now I shrink the vector on order to optimize
  // the memory consnmption.
  tripletList.shrink_to_fit();

  Eigen::SparseMatrix<geom::real> curA(dim_, dim_);

  curA.setFromTriplets(tripletList.begin(), tripletList.end());

  // I remove numerical zeros, I hope it works well
  curA.prune(A_.coeff(0,0));

  if(sym_ == symExpr)
    A_ += curA;
  else if(sym_ == false && symExpr == true)
    A_ += curA.selfadjointView<Eigen::Upper>();
  else
  {
    std::cerr << "The symmetry has been broken" << std::endl;
    sym_ = false;
    A_ = A_.selfadjointView<Eigen::Upper>();
    A_ += curA;
  }
}

template <typename T>
void Assembler::integrateVolRhs(const ExprWrapper<T>& expr)
{
  const T& exprDerived(expr);

  unsigned elemNo = 0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = elemNo * Vh_.getDofNo();

    for(unsigned i = 0; i < Vh_.getDofNo(); i++)
      for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
          b_(i + indexOffset) += exprDerived(*it, i, t, q) *
                                 it->getWeight(q) *
                                 it->getAbsDetJac(t);
    elemNo++;
  }

}

template <typename T>
void Assembler::integrateFacesExtRhs(const ExprWrapper<T>& expr, unsigned BClabel)
{
  const T& exprDerived(expr);

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == BClabel)
    {
      unsigned indexOffset = it->getElem() * Vh_.getDofNo();
      for(unsigned i = 0; i < Vh_.getDofNo(); i++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
          b_(i + indexOffset) += exprDerived(*it, i, q) *
                                 it->getWeight(q) *
                                 it->getAreaDoubled();
    }
}

} // namespace dgfem

#endif // _ASSEMBLER_HPP_
