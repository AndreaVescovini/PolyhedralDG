/*!
    @file   Problem.hpp
    @author Andrea Vescovini
    @brief  Class for solving a differential problem
*/

#ifndef _PROBLEM_HPP_
#define _PROBLEM_HPP_

#include "ExprWrapper.hpp"
#include "FeSpace.hpp"
#include "PolyDG.hpp"
#include "Watch.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace PolyDG
{

/*!
    @brief  Class for solving a differential problem

    This class is the key one in the solution of a differential problem through
    discontinuous finite elements.@n
    The key methods are integrateVol, integrateFacesExt, integrateFacesInt that
    allow to integrate the expression of a bilinear form over the volume, external
    faces and integral faces of the elements of the mesh. The method receive
    directly the expression related to the bilinear form of the problem and
    compute the integrals that are needed to assembleme the matrix of the linear
    system through the tecnique of expression templates. Once all the integrate
    functions have been called, the method finalizeMatrix() is needed in order
    to effecively fill the matrix.@n
    Then in an analogous way through integrateVolRhs and integrateFacesExtRhs
    the integrals of the linear functional can be performed, in order to allow
    the assembling of the rhs of the linear system.

    After that the linear system is ready to be solved using a direct method
    (LU decomposition or Chlolesky decomposition) or an iterative one (conjugate
    gradient or BiCGSTAB), depending also on the symmetry or non-symmetry of the
    problem. If you want you can get the matrix and the rhs and solve the system
    with your own method.

    Finally the errors can be computed and the solution exported into a VTK file,
    in order to be read with a visualization software (i.e. Paraview).@n
    If you want to use the same Problem to solve another problem you must call
    the method clearMatrix() and clearRhs() in order to clean the containers.
*/

class Problem
{
public:
  //! Constructor
  explicit Problem(const FeSpace& Vh);

  //! Copy constructor
  Problem(const Problem&) = default;

  //! Move constructor
  Problem(Problem&&) = default;

  /*!
      @brief Integrate the expression of a bilinear form over the volume

      This function integrates the expression of a bilinear form over the volume
      of the elements of the mesh, in order to assemble the matrix of the linear
      system. The whole expressione is assembled through the tecnique of
      expression templates. If the bilinear form is symmetric it is conveninet
      to specify sym = true in order to compute and store only half of the
      integrals (those related to the upper triangular part of the matrix).

      @param expr Expression of a bilinear form.
      @param sym  @c true if the form is symmetric, @c false if it is not.
  */
  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool sym = false);

  /*!
      @brief Integrate an expression of a bilinear form over the external faces

      This function integrates the expression of a bilinear form over the external
      faces of the mesh the match the given BCLabelType, in order to assemble the
      matrix of the linear system. The whole expressione is assembled through
      the tecnique of expression templates. If the bilinear form is symmetric it
      is conveninet to specify sym = true in order to compute and store only
      half of the integrals (those related to the upper triangular part of the
      matrix).

      @param expr     Expression of a bilinear form.
      @param bcLabels Vector of BCLabelType used to specify over which external
                      faces the integration has to be done.
      @param sym      @c true if the form is symmetric, @c false if it is not.
  */
  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels,
                         bool sym = false);

  /*!
      @brief Integrate an expression of a bilinear form over the internal faces

      This function integrates the expression of a blinear form over the internal
      faces of the mesh, in order to assemble the matrix of the linear system.
      The whole expressione is assembled through the tecnique of expression
      templates. If the bilinear form is symmetric it is conveninet to specify
      sym = true in order to compute and store only half of the integrals
      (those related to the upper triangular part of the matrix).

      @param expr Expression of a bilinear form.
      @param sym  @c true if the form is symmetric, @c false if it is not.
  */
  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool sym = false);

  /*!
      @brief Integrate an expression of a linear form over the volume

      This function integrates the expression of a linear form over the volume
      of the elements of the mesh the match the given BCLabelType, in order to assemble
      the rhs of the linear system. The whole expressione is assembled through
      the tecnique of expression templates.

      @param expr Expression of a linear form.
  */
  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);

  /*!
      @brief Integrate an expression of a linear form over the external faces

      This function integrates the expression of a linear form over the external
      faces of the mesh, in order to assmble the rhs of the linear system. The
      whole expressione is assembled through the tecnique of expression templates.

      @param expr     Expression of a linear form.
      @param bcLabels Vector of BCLabelType used to specify over which external
                      faces the integration has to be done.
  */
  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels);

  /*!
      @brief Solve the linear system with a sparse LU decomposition

      Thi function solves the linear system with a direct sparse LU decomposition,
      using the solver implemented in the library Eigen (https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html).
      The solver works both for symmetric and non symmetric matrices.
      If the decomposition fails, a @c std::runtime_error exception is thrown.
      If the solver does not succeed a warning message is printed.

      @return @c true if the solver succeeds, @c false if it does not.
  */
  bool solveLU();

  /*!
      @brief Solve the linear system with a sparse Chlolesky decomposition

      Thi function solves the linear system with a direct sparse Chlolesky
      decomposition, using the solver implemented in the library Eigen
      (https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLLT.html). The
      solver works only for symmetric matrices, it throws a @c std::domain_error
      exception if called on a non symmetric matrix. If the decomposition fails,
      a @c std::runtime_error exception is thrown. If the solver does not succeed
      a warning message is printed.

      @return @c true if the solver succeeds, @c false if it does not.
  */
  bool solveCholesky();

  /*!
      @brief Solve the linear system with the conjugate gradient method

      This function solves the linear system iteratively using the conjugate
      gradient method, using the solver implemented in the library Eigen
      (https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html).The
      solver works only for symmetric matrices, it throws a @c std::domain_error
      exception if called on a non symmetric matrix. It requires and initial guess
      and you can specify the maximum number of iteation and the tolerance for
      the stopping criterion. If it does not converge in the assigned maximum
      number of iterarions a warning message is printed.

      @param x0      Initial guess.
      @param iterMax Maximum number if iteration, if not specified it is 10000.
      @param tol     Tolerance for the stopping criterio, if not specified it is
                     the machine epsilon.
      @return @c true if the solver succeeds, @c false if it does not.

  */
  bool solveCG(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
               Real tol = Eigen::NumTraits<Real>::epsilon());

  /*!
      @brief Solve the linear system with the bi-conjugate gradient stabilized gradient method

      This function solves the linear system iteratively using the bi-conjugate
      gradient stabilized method, using the solver implemented in the library Eigen
      (https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html).The solver
      works for symmetric and non symmetric matrices, but it convenient only for
      non symmetric matrices. It requires and initial guess and you can specify
      the maximum number of iteation and the tolerance for the stopping criterion.
      If it does not converge in the assigned maximum number of iterarions a
      warning message is printed.

      @param x0      Initial guess.
      @param iterMax Maximum number if iteration, if not specified it is 10000.
      @param tol     Tolerance for the stopping criterio, if not specified it is
                     the machine epsilon.
      @return @c true if the solver succeeds, @c false if it does not.
  */
  bool solveBiCGSTAB(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
                     Real tol = Eigen::NumTraits<Real>::epsilon());

  /*!
      @brief Compute the L-2 norm of the error

      This function, given the exact solution uex, computes \f$ || u_h - u_{ex} ||_{L^2(\mathcal{T})} \f$.

      @param uex Exact solution, it is a function that takes a @c Eigen::Vector3d
                 and gives a PolyDG::Real.
  */
  Real computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const;

  /*!
      @brief Compute the H1-seminorm of the error

      This function, given the gradient of the exact solution uexGrad, computes
      \f$ || \nabla u_h - \nabla u_{ex} ||_{L^2(\mathcal{T})} \f$.

      @param uexGrad Gradient of the exact solution, it is a function that takes a
                     @c Eigen::Vector3d and gives a @c Eigen::Vector3d.
  */
  Real computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const;

  /*!
      @brief Export the solution

      This function exports the solution into a VTK unstructured grid file with
      XML format. It can be read with a visualization software (e.g. Paraview).

      @param fileName  Name of the file to be saved (the extension should be .vtu).
      @param precision Precision to be used for floating points numbers.
  */
  void exportSolutionVTK(const std::string& fileName, unsigned precision = 8) const;

  /*!
      @brief Export the solution

      This function exports the vector u (supposed to be the solution obtained
      with an external solver) into a VTK unstructured grid file with XML
      format. It can be read with a visualization software (e.g. Paraview).

      @param u         Vector containing the solution of the problem.
      @param fileName  Name of the file to be saved (the extension should be .vtu).
      @param precision Precision to be used for floating points numbers.
  */
  void exportSolutionVTK(const Eigen::VectorXd& u, const std::string& fileName, unsigned precision = 8) const;

  /*!
      @brief Get the symmetry

      This function tells if globally the variational form and so the matrix is
      symmetric or not.

      @return @c true if the system is symmetric, @c false if it is not.
  */
  bool isSymmetric() const;

  //! Get the sparse matrix of the linear system
  inline const Eigen::SparseMatrix<Real>& getMatrix() const;

  //! Get the rhs of the linear system
  inline const Eigen::VectorXd& getRhs() const;

  //! Get the solution of the linear system
  inline const Eigen::VectorXd& getSolution() const;

  //! Get the dimension of the linear system
  inline unsigned getDim() const;

  /*!
      @brief Clear the matrix of the linear system

      This function sets to zero the matrix of the linear system. It has to be
      called if you want to integrate a new problem.
  */
  void clearMatrix();

  /*!
      @brief Clear the rhs of the linear system

      This function sets to zero the rhs of the linear system. It has to be
      called if you want to integrate a new problem.
  */
  void clearRhs();

  /*!
      @brief Assemble the matrix of the linear system.

      This function assembles the matrix of the linear system. It has to be
      after the integration methods and before the solve methods.
  */
  void finalizeMatrix();

  //! Print information about the problem
  void printInfo(std::ostream& out = std::cout) const;

  //! Destructor
  virtual ~Problem() = default;

private:
  //! Alias for a triplet (row, column, value)
  using triplet = Eigen::Triplet<double>;

  //! FeSpace used for the solution of the problem
  const FeSpace& Vh_;

  //! Dimension of the linear system
  const unsigned dim_;

  //! Matrix of the linear system
  Eigen::SparseMatrix<Real> A_;

  //! Rhs of the linear system
  Eigen::VectorXd b_;

  //! Solution of the linear system
  Eigen::VectorXd u_;

  //! Containers used to store the triplets computed in the integration
  std::vector<std::vector<triplet>> triplets_;

  //! Vector of bools referring to the symmetry of the integrated forms
  std::vector<bool> sym_;

  /*!
      @brief Evaluate the solution
      @param u  The vector containing the solution.
      @param x  /f$ x /f$ coordinate in the element el at which the solution has
                to be evaluated.
      @param y  /f$ y /f$ coordinate in the element el at which the solution has
                to be evaluated.
      @param z  /f$ z /f$ coordinate in the element el at which the solution has
                to be evaluated.
      @param el FeElement in which the solution has to be evaluated.
  */
  Real evalSolution(const Eigen::VectorXd& u, Real x, Real y, Real z, const FeElement& el) const;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename T>
void Problem::integrateVol(const ExprWrapper<T>& expr, bool sym)
{
  #ifdef VERBOSITY
    std::cout << "Integrating the bilinear form over volumes..............";
    Utilities::Watch ch;
    ch.start();
  #endif

  // I exploit the conversion to derived
  const T& exprDerived(expr);

  triplets_.emplace_back();
  sym_.push_back(sym);

  // If the variational form is symmetric I store only half elements
  if(sym == true)
    triplets_.back().reserve(dim_ * (Vh_.getDof() + 1) / 2);
  else
    triplets_.back().reserve(dim_ * Vh_.getDof());

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    const unsigned indexOffset = it->getElem().getId() * Vh_.getDof();

    for(unsigned j = 0; j < Vh_.getDof(); j++)
      for(unsigned i = 0; i < (sym == true ? j + 1 : Vh_.getDof()); i++)
      {
        Real sum = 0.0;
        for(SizeType t = 0; t < it->getTetrahedraNo(); t++)
          for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
            sum += exprDerived(*it, i, j, t, p) * it->getWeight(p) * it->getAbsDetJac(t);

        triplets_.back().emplace_back(i + indexOffset, j + indexOffset, sum);
      }
  }

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
}

template <typename T>
void Problem::integrateFacesExt(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels,
                                bool sym)
{
  #ifdef VERBOSITY
    std::cout << "Integrating the bilinear form over external faces (";
    for(SizeType i = 0; i < bcLabels.size() - 1; i++)
      std::cout << bcLabels[i] << ", ";

    std::cout << bcLabels.back() << ")...";
    Utilities::Watch ch;
    ch.start();
  #endif

  const T& exprDerived(expr);

  triplets_.emplace_back();
  sym_.push_back(sym);

  // If the variational form is symmetric I store only half elements,
  // I overestimate considering all the external faces with the same type of
  // boundary conditions
  if(sym == true)
    triplets_.back().reserve(Vh_.getFeFacesExtNo() * Vh_.getDof() * (Vh_.getDof() + 1) / 2);
  else
    triplets_.back().reserve(Vh_.getFeFacesExtNo() * Vh_.getDof() * Vh_.getDof());

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(std::find(bcLabels.cbegin(), bcLabels.cend(), it->getBClabel()) != bcLabels.cend())
    {
      const unsigned indexOffset = it->getElemIn() * Vh_.getDof();

      for(unsigned j = 0; j < Vh_.getDof(); j++)
        for(unsigned i = 0; i < (sym == true ? j + 1 : Vh_.getDof()); i++)
        {
          Real sum = 0.0;

          for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
            sum += exprDerived(*it, i, j, p) * it->getWeight(p) * it->getAreaDoubled();

          triplets_.back().emplace_back(i + indexOffset, j + indexOffset, sum);
        }
    }

  // I made and overestimate so now I shrink the vector in order to optimize
  // the memory consnmption.
  triplets_.back().shrink_to_fit();

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
}

template <typename T>
void Problem::integrateFacesInt(const ExprWrapper<T>& expr, bool sym)
{
  #ifdef VERBOSITY
    std::cout << "Integrating the bilinear form over internal faces.......";
    Utilities::Watch ch;
    ch.start();
  #endif

  const T& exprDerived(expr);

  triplets_.emplace_back();
  sym_.push_back(sym);

  // If the variational form is symmetric I store only half elements
  if(sym == true)
    triplets_.back().reserve(Vh_.getFeFacesIntNo() * Vh_.getDof() * (2 * Vh_.getDof() + 1));
  else
    triplets_.back().reserve(Vh_.getFeFacesIntNo() * Vh_.getDof() * Vh_.getDof() * 4);

  const std::array<SideType, 2> sides = {{Out, In}};

  for(auto it = Vh_.feFacesIntCbegin(); it != Vh_.feFacesIntCend(); it++)
  {
    const std::array<unsigned, 2> indexOffset = {{ it->getElemIn() * Vh_.getDof(),
                                                   it->getElemOut() * Vh_.getDof() }};

    for(unsigned sj = 0; sj < 2; sj++)
      for(unsigned si = 0; si < (sym == true ? sj + 1 : 2); si++)
        for(unsigned j = 0; j < Vh_.getDof(); j++)
          for(unsigned i = 0; i < (sym == true && sides[si] == sides[sj] ? j + 1 : Vh_.getDof()); i++)
          {
            Real sum = 0.0;

            for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
              sum += exprDerived(*it, i, j, sides[si], sides[sj], p) * it->getWeight(p) * it->getAreaDoubled();

            triplets_.back().emplace_back(i + indexOffset[si], j + indexOffset[sj], sum);
          }
  }

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
}

template <typename T>
void Problem::integrateVolRhs(const ExprWrapper<T>& expr)
{
  #ifdef VERBOSITY
    std::cout << "Integrating the rhs over volumes........................";
    Utilities::Watch ch;
    ch.start();
  #endif

  const T& exprDerived(expr);

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    const unsigned indexOffset = it->getElem().getId() * Vh_.getDof();

    for(unsigned i = 0; i < Vh_.getDof(); i++)
      for(SizeType t = 0; t < it->getTetrahedraNo(); t++)
        for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
          b_(i + indexOffset) += exprDerived(*it, i, t, p) *
                                 it->getWeight(p) *
                                 it->getAbsDetJac(t);
  }

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
}

template <typename T>
void Problem::integrateFacesExtRhs(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels)
{
  #ifdef VERBOSITY
    std::cout << "Integrating the rhs over external faces (";
    for(SizeType i = 0; i < bcLabels.size() - 1; i++)
      std::cout << bcLabels[i] << ", ";

    std::cout << bcLabels.back() << ").............";
    Utilities::Watch ch;
    ch.start();
  #endif

  const T& exprDerived(expr);

  for(auto it = Vh_.feFacesExtCbegin(); it != Vh_.feFacesExtCend(); it++)
    if(std::find(bcLabels.cbegin(), bcLabels.cend(), it->getBClabel()) != bcLabels.cend())
    {
      unsigned indexOffset = it->getElemIn() * Vh_.getDof();
      for(unsigned i = 0; i < Vh_.getDof(); i++)
        for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
        {
          b_(i + indexOffset) += exprDerived(*it, i, p) *
                                 it->getWeight(p) *
                                 it->getAreaDoubled();
        }
    }

  #ifdef VERBOSITY
    ch.stop();
    std::cout << "Done!   " << ch << std::endl;
  #endif
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

} // namespace PolyDG

#endif // _PROBLEM_HPP_
