#include "Assembler.hpp"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

namespace dgfem
{

Assembler::Assembler(const FeSpace& Vh, bool sym)
  : Vh_{Vh}, dim_{Vh.getFeElementsNo() * Vh.getDofNo()}, A_{dim_, dim_},
    b_{Eigen::VectorXd::Zero(dim_)}, u_{Eigen::VectorXd::Zero(dim_)}, sym_{sym} {}

void Assembler::solveLU()
{
  Eigen::SparseLU<Eigen::SparseMatrix<geom::real>> solver;

  if(sym_ == true)
  {
    solver.isSymmetric(true);
    solver.compute(A_.selfadjointView<Eigen::Upper>());
  }
  else
    solver.compute(A_);

  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue\n" << solver.lastErrorMessage() << std::endl;
    return;
  }

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue\n" << solver.lastErrorMessage() << std::endl;
    return;
  }
}

void Assembler::solveChol()
{
  if(sym_ == false)
  {
    std::cerr << "solveChol() requires a symmetric matrix" << std::endl;
    return;
  }

  Eigen::SimplicialLLT<Eigen::SparseMatrix<geom::real>, Eigen::Upper> solver;
  // Eigen::SimplicialLDLT<Eigen::SparseMatrix<geom::real>, Eigen::Upper> solver;

  A_.makeCompressed();
  solver.compute(A_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue " << solver.info() << std::endl;
    return;
  }

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue" << std::endl;
    return;
  }
}

void Assembler::solveCG(const Eigen::VectorXd& x0, unsigned iterMax, geom::real tol)
{
  if(sym_ == false)
  {
    std::cerr << "solveCG() requires a symmetric matrix" << std::endl;
    return;
  }

  Eigen::ConjugateGradient<Eigen::SparseMatrix<geom::real>, Eigen::Upper> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<geom::real>, Eigen::Lower|Upper> solver;

  solver.setMaxIterations(iterMax);
  solver.setTolerance(tol);

  solver.compute(A_);

  u_ = solver.solveWithGuess(b_, x0);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Not converged" << std::endl;
    return;
  }

  // if(verbosity)
  std::cout << "Converged with " << solver.iterations() << " iterations.\n";
  std::cout << "Estimated error "<< solver.error() << std::endl;
}

void Assembler::isSymmetric(bool sym)
{
  sym_ = sym;
}


const Eigen::SparseMatrix<geom::real> Assembler::getMatrix() const
{
  return A_;
}

const Eigen::VectorXd Assembler::getRhs() const
{
  return b_;
}

const Eigen::VectorXd Assembler::getSolution() const
{
  return u_;
}

unsigned Assembler::getDim() const
{
  return dim_;
}

void Assembler::clearMatrix()
{
  A_.setZero();
}

void Assembler::clearRhs()
{
  b_ = Eigen::VectorXd::Zero(dim_);
}

void Assembler::printMatrix(std::ostream& out) const
{
  out << A_ << std::endl;
}

void Assembler::printMatrixSym(std::ostream& out) const
{
  out << A_.selfadjointView<Eigen::Upper>() << std::endl;
}

void Assembler::printRhs(std::ostream& out) const
{
  out << b_ << std::endl;
}


} // namespace dgfem
