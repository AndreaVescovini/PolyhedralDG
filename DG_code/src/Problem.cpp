#include "Problem.hpp"
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <cmath>

namespace dgfem
{

Problem::Problem(const FeSpace& Vh, bool sym)
  : Vh_{Vh}, dim_{Vh.getFeElementsNo() * Vh.getDofNo()}, A_{dim_, dim_},
    b_{Eigen::VectorXd::Zero(dim_)}, u_{Eigen::VectorXd::Zero(dim_)}, sym_{sym} {}

void Problem::solveLU()
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

void Problem::solveChol()
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

void Problem::solveCG(const Eigen::VectorXd& x0, unsigned iterMax, geom::real tol)
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

geom::real Problem::computeErrorL2(const std::function<geom::real (Eigen::Vector3d)>& uex) const
{
  geom::real errSquared = 0.0;

  unsigned elemNo = 0;
  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = elemNo * Vh_.getDofNo();
    for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
      for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
      {
        geom::real uh = 0.0;

        // Evaluation of the fem function at the quadrature node
        for(unsigned i = 0; i < Vh_.getDofNo(); i++)
          uh += u_(i + indexOffset) * it->getPhi(t, q, i);

        geom::real difference = uh - uex(it->getQuadPoint(t, q));
        errSquared += difference * difference * it->getWeight(q) * it->getAbsDetJac(t);
      }

    elemNo++;
  }

  return std::sqrt(errSquared);
}

geom::real Problem::computeErrorH10(const std::function<Eigen::Vector3d (Eigen::Vector3d)>& uexGrad) const
{
  geom::real errSquared = 0.0;

  unsigned elemNo = 0;
  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = elemNo * Vh_.getDofNo();
    for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
      for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
      {
        Eigen::Vector3d uhGrad = Eigen::Vector3d::Zero();

        // Evaluation of the fem function at the quadrature node
        for(unsigned i = 0; i < Vh_.getDofNo(); i++)
          uhGrad += u_(i + indexOffset) * it->getPhiDer(t, q, i);

        Eigen::Vector3d difference = uhGrad - uexGrad(it->getQuadPoint(t, q));
        errSquared += difference.squaredNorm() * it->getWeight(q) * it->getAbsDetJac(t);
      }

    elemNo++;
  }

  return std::sqrt(errSquared);
}

void Problem::isSymmetric(bool sym)
{
  sym_ = sym;
}

bool Problem::getSymmetry() const
{
  return sym_;
}

const Eigen::SparseMatrix<geom::real> Problem::getMatrix() const
{
  return A_;
}

const Eigen::VectorXd Problem::getRhs() const
{
  return b_;
}

const Eigen::VectorXd Problem::getSolution() const
{
  return u_;
}

unsigned Problem::getDim() const
{
  return dim_;
}

void Problem::clearMatrix()
{
  A_.setZero();
}

void Problem::clearRhs()
{
  b_ = Eigen::VectorXd::Zero(dim_);
}

} // namespace dgfem
