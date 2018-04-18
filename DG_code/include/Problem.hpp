#ifndef _PROBLEM_HPP_
#define _PROBLEM_HPP_

#include "PolyDG.hpp"
#include "FeSpace.hpp"
#include "ExprWrapper.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <string>

namespace PolyDG
{

class Problem
{
public:
  explicit Problem(const FeSpace& Vh, bool sym = false);

  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool symExpr = false);

  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool symExpr = false);

  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr, BCtype bcLabel, bool symExpr = false);

  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);

  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr, BCtype bcLabel);

  void solveLU();
  void solveChol();
  void solveCG(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
               Real tol = Eigen::NumTraits<Real>::epsilon());

  Real computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const;
  Real computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const;

  void exportSolutionVTK(const std::string& fileName) const;

  inline void isSymmetric(bool sym);
  inline bool getSymmetry() const;

  inline const Eigen::SparseMatrix<Real> getMatrix() const;
  inline const Eigen::VectorXd getRhs() const;
  inline const Eigen::VectorXd getSolution() const;
  inline unsigned getDim() const;

  inline void clearMatrix();
  inline void clearRhs();

  virtual ~Problem() = default;

private:
  const FeSpace& Vh_;
  unsigned dim_;
  Eigen::SparseMatrix<Real> A_;
  // std::vector<triplet> tripletList
  Eigen::VectorXd b_;
  Eigen::VectorXd u_;

  bool sym_;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename T>
void Problem::integrateVol(const ExprWrapper<T>& expr, bool symExpr)
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
        Real sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        tripletList.emplace_back(i + indexOffset, j + indexOffset, sum);
      }

    elemNo++;
  }

  Eigen::SparseMatrix<Real> curA(dim_, dim_);

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
void Problem::integrateFacesInt(const ExprWrapper<T>& expr, bool symExpr)
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

            Real sum = 0.0;

            for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
              sum += exprDerived(*it, i, j, side1, side2, q) * it->getWeight(q) * it->getAreaDoubled();

            unsigned rowIndex = (symExpr == true ? std::min(i + indexOffset[side1], j + indexOffset[side2]) : i + indexOffset[side1] );
            unsigned colIndex = (symExpr == true ? std::max(i + indexOffset[side1], j + indexOffset[side2]) : j + indexOffset[side2] );

            tripletList.emplace_back(rowIndex, colIndex, sum);
          }
  }

  Eigen::SparseMatrix<Real> curA(dim_, dim_);

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
void Problem::integrateFacesExt(const ExprWrapper<T>& expr, BCtype bcLabel, bool symExpr)
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
    if(it->getBClabel() == bcLabel)
    {
      unsigned indexOffset = it->getElem() * Vh_.getDofNo();

      for(unsigned j = 0; j < Vh_.getDofNo(); j++)
        for(unsigned i = 0; i < (symExpr == true ? j + 1 : Vh_.getDofNo()); i++)
        {
          Real sum = 0.0;

          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, q) * it->getWeight(q) * it->getAreaDoubled();

          tripletList.emplace_back(i + indexOffset, j + indexOffset, sum);
        }
    }

  // I made and overestimate so now I shrink the vector on order to optimize
  // the memory consnmption.
  tripletList.shrink_to_fit();

  Eigen::SparseMatrix<Real> curA(dim_, dim_);

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
void Problem::integrateVolRhs(const ExprWrapper<T>& expr)
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
void Problem::integrateFacesExtRhs(const ExprWrapper<T>& expr, BCtype bcLabel)
{
  const T& exprDerived(expr);

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == bcLabel)
    {
      unsigned indexOffset = it->getElem() * Vh_.getDofNo();
      for(unsigned i = 0; i < Vh_.getDofNo(); i++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
        {
          b_(i + indexOffset) += exprDerived(*it, i, q) *
                                 it->getWeight(q) *
                                 it->getAreaDoubled();
        }
    }
}

inline bool Problem::getSymmetry() const
{
  return sym_;
}

inline void Problem::isSymmetric(bool sym)
{
  sym_ = sym;
}

inline const Eigen::SparseMatrix<Real> Problem::getMatrix() const
{
  return A_;
}

inline const Eigen::VectorXd Problem::getRhs() const
{
  return b_;
}

inline const Eigen::VectorXd Problem::getSolution() const
{
  return u_;
}

inline unsigned Problem::getDim() const
{
  return dim_;
}

inline void Problem::clearMatrix()
{
  A_.setZero();
}

inline void Problem::clearRhs()
{
  b_ = Eigen::VectorXd::Zero(dim_);
}

} // namespace PolyDG

#endif // _PROBLEM_HPP_
