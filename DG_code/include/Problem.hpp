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
#include <array>

namespace PolyDG
{

class Problem
{
public:

  // Constructor that requires a FeSpace used for the computation of the solution
  // and for the test functions.
  explicit Problem(const FeSpace& Vh);

  // Function that integrates a blinear form over the volume of the
  // elements of the mesh. If sym == true only the upper triangular part is
  // computed.
  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool sym = false);

  // Function that integrates a bilinear form over the internal faces
  // of the elements of the mesh. If sym == true only the upper triangular part
  // is computed.
  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool sym = false);

  // Function that integrates a bilinear over the external faces
  // of the elements of the mesh. If sym == true only the upper triangular part
  // is computed.
  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr, BCtype bcLabel, bool sym = false);

  // Function that integrates a linear form over the volume of the elements of
  // the mesh.
  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);

  // Function that integrates a linear form over the external faces of the
  // elements of the mesh.
  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr, BCtype bcLabel);

  // Function that solves the linear system with a direct LU decomposition.
  void solveLU();

  // Function that solves the linear system with a direct Chlolesky decomposition.
  // It works only with a symmetric formulation.
  void solveChol();

  // Function that solves the linear system iteratively using the Conjugate
  // Gradient method. It works only with a symmetric formulation. it requires
  // an initial guess and maximum number of iterations and the tolerance can be
  // specified.
  void solveCG(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
               Real tol = Eigen::NumTraits<Real>::epsilon());

  // Function that, given the exact solution uex, computes the L2-norm of the
  // error.
  Real computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const;

  // Function that, given the gradient of the exact solution uexGrad, computes
  // the H1-seminorm of the error.
  Real computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const;

  // Function that exports the solution in a xml file fileName.vtu that can be
  // read with a visualization software (e.g. Paraview).
  void exportSolutionVTK(const std::string& fileName) const;

  // Function that returns true if the inserted form is symmetric or false
  // if it is not.
  inline bool isSymmetric() const;

  // Function that returns a const reference to the matrix of the linear system.
  inline const Eigen::SparseMatrix<Real>& getMatrix() const;

  // Functions that returns a const reference to the vector of the rhs of the
  // linear system.
  inline const Eigen::VectorXd& getRhs() const;

  // Functions that returns a const reference to the vector of the solution of the
  // linear system.
  inline const Eigen::VectorXd& getSolution() const;

  // Function that return the dimension of the linear system.
  inline unsigned getDim() const;

  // Function that sets to zero the matrix of the linear system.
  inline void clearMatrix();

  // Function that sets to zero the vector of the rhs of the linear system.
  inline void clearRhs();

  // Function that assembles the matrix of the linear system. This must be called
  // after the integrations.
  void finalizeMatrix();

  virtual ~Problem() = default;

private:
  using triplet = Eigen::Triplet<double>;

  // FeSpace used for the solution of the problem and for the test functions.
  const FeSpace& Vh_;

  // Dimension of the linear system.
  unsigned dim_;

  // Matrix, rhs and solution of the linear system.
  Eigen::SparseMatrix<Real> A_;
  Eigen::VectorXd b_;
  Eigen::VectorXd u_;

  // Array of vectors of triplet used to store the values computed in the
  // integration.
  std::array<std::vector<triplet>, 3> triplets_;

  // Array of bools referring to the symmetriy of the bilinear forms respectively
  // over volumes, internal faces and external faces.
  std::array<bool, 3> sym_;

  // Function that evalues the solution of the problem at the point x, y, z
  // belonging to the FeElement el.
  Real evalSolution(Real x, Real y, Real z, const FeElement& el) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename T>
void Problem::integrateVol(const ExprWrapper<T>& expr, bool sym)
{
  // I exploit the conversion to derived
  const T& exprDerived(expr);

  sym_[0] = sym;

  // If the variational form is symmetric I store only half elements
  if(sym == true)
    triplets_[0].reserve(dim_ * (Vh_.getDofNo() + 1) / 2);
  else
    triplets_[0].reserve(dim_ * Vh_.getDofNo());

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = it->getElem().getId() * Vh_.getDofNo();

    for(unsigned j = 0; j < Vh_.getDofNo(); j++)
      for(unsigned i = 0; i < (sym == true ? j + 1 : Vh_.getDofNo()); i++)
      {
        Real sum = 0.0;
        for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, t, q) * it->getWeight(q) * it->getAbsDetJac(t);

        triplets_[0].emplace_back(i + indexOffset, j + indexOffset, sum);
      }
  }
}

template <typename T>
void Problem::integrateFacesInt(const ExprWrapper<T>& expr, bool sym)
{
  const T& exprDerived(expr);

  sym_[1] = sym;

  // If the variational form is symmetric I store only half elements
  if(sym == true)
    triplets_[1].reserve(Vh_.getFeFacesIntNo() * Vh_.getDofNo() * (2 * Vh_.getDofNo() + 1));
  else
    triplets_[1].reserve(Vh_.getFeFacesIntNo() * Vh_.getDofNo() * Vh_.getDofNo() * 4);

  for(auto it = Vh_.feFacesIntCbegin(); it != Vh_.feFacesIntCend(); it++)
  {
    std::array<unsigned, 2> indexOffset = { it->getElem(0) * Vh_.getDofNo(),
                                            it->getElem(1) * Vh_.getDofNo() };

    for(int side2 = 0; side2 < 2; side2++)
      for(int side1 = 0; side1 < (sym == true ? side2 + 1 : 2); side1++)
        for(unsigned j = 0; j < Vh_.getDofNo(); j++)
          for(unsigned i = 0; i < (sym == true && side1 == side2 ? j + 1 : Vh_.getDofNo()); i++)
          {
            Real sum = 0.0;

            for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
              sum += exprDerived(*it, i, j, side1, side2, q) * it->getWeight(q) * it->getAreaDoubled();

            triplets_[1].emplace_back(i + indexOffset[side1], j + indexOffset[side2], sum);
          }
  }
}

template <typename T>
void Problem::integrateFacesExt(const ExprWrapper<T>& expr, BCtype bcLabel, bool sym)
{
  const T& exprDerived(expr);

  sym_[2] = sym;

  // If the variational form is symmetric I store only half elements,
  // I overestimate considering all the external faces with the same type of
  // boundary conditions
  if(sym == true)
    triplets_[2].reserve(Vh_.getFeFacesExtNo() * Vh_.getDofNo() * (Vh_.getDofNo() + 1) / 2);
  else
    triplets_[2].reserve(Vh_.getFeFacesExtNo() * Vh_.getDofNo() * Vh_.getDofNo());

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(it->getBClabel() == bcLabel)
    {
      unsigned indexOffset = it->getElem() * Vh_.getDofNo();

      for(unsigned j = 0; j < Vh_.getDofNo(); j++)
        for(unsigned i = 0; i < (sym == true ? j + 1 : Vh_.getDofNo()); i++)
        {
          Real sum = 0.0;

          for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
            sum += exprDerived(*it, i, j, q) * it->getWeight(q) * it->getAreaDoubled();

          triplets_[2].emplace_back(i + indexOffset, j + indexOffset, sum);
        }
    }

  // I made and overestimate so now I shrink the vector on order to optimize
  // the memory consnmption.
  triplets_[2].shrink_to_fit();
}

template <typename T>
void Problem::integrateVolRhs(const ExprWrapper<T>& expr)
{
  const T& exprDerived(expr);

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = it->getElem().getId() * Vh_.getDofNo();

    for(unsigned i = 0; i < Vh_.getDofNo(); i++)
      for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
        for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
          b_(i + indexOffset) += exprDerived(*it, i, t, q) *
                                 it->getWeight(q) *
                                 it->getAbsDetJac(t);
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

inline bool Problem::isSymmetric() const
{
  return (sym_[0] || sym_[1] || sym_[2]);
}

inline const Eigen::SparseMatrix<Real>& Problem::getMatrix() const
{
  return A_;
}

inline const Eigen::VectorXd& Problem::getRhs() const
{
  return b_;
}

inline const Eigen::VectorXd& Problem::getSolution() const
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
