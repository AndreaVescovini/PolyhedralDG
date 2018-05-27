/*!
    @file   test_solver.cpp
    @author Andrea Vescovini
    @brief  Test for the solvers of the linear system
*/

// Test for the solvers of the linear system.
//
// 1)  -laplacian(u) = f   in omega
//                 u = gd   in delta_omega
//
// 2)  -laplacian(u) + u = f   in omega
//                     u = gd  in delta_omega

#include "ExprOperators.hpp"
#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "PolyDG.hpp"
#include "Problem.hpp"
#include "Utilities.hpp"
#include "Watch.hpp"

#include <Eigen/Core>
#include "GetPot.hpp"

#include <cmath>
#include <exception>
#include <string>
#include <vector>

/*!
    Two problems are solved with the available methods:

    1)  \f$ - \Delta u = f  \quad \text{in} \quad \Omega\\
                     u = g_d  \quad \text{on} \quad \partial \Omega \f$ @n
      using a symmetric formulation and:

    2)  \f$ - \Delta u + u = f  \quad \text{in} \quad \Omega\\
                     u = g_d  \quad \text{on} \quad \partial \Omega \f$ @n
      using a non-symmetric formulation.
*/

int main(int argc, char* argv[])
{
  using Utilities::pow;

  // Exact solution and source term
  auto uex = [](const Eigen::Vector3d& x) { return std::exp(x(0) * x(1) * x(2)); };
  auto source = [&uex](const Eigen::Vector3d& x) { return -uex(x) * (pow(x(0) * x(1), 2) +
                                                                     pow(x(1) * x(2), 2) +
                                                                     pow(x(0) * x(2), 2));};
  auto uexGrad = [&uex](const Eigen::Vector3d& x) -> Eigen::Vector3d { return uex(x)*Eigen::Vector3d(x(1) * x(2),
                                                                                                     x(0) * x(2),
                                                                                                     x(0) * x(1)); };

 GetPot comLine(argc, argv);
 const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
 GetPot fileData(fileName.c_str());

 const std::string meshFile = fileData("dir", "../../meshes") + "/cube_str1296p.mesh";

  // Mesh reading
  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(meshFile, reader);
  Th.printInfo();

  // FeSpace Creation
  unsigned r = 2;
  PolyDG::FeSpace Vh(Th, r);

  // Operators
  PolyDG::PhiI            v;
  PolyDG::GradPhiJ        uGrad;
  PolyDG::GradPhiI        vGrad;
  PolyDG::JumpPhiJ        uJump;
  PolyDG::JumpPhiI        vJump;
  PolyDG::AverGradPhiJ    uGradAver;
  PolyDG::AverGradPhiI    vGradAver;
  PolyDG::PenaltyScaling  gamma(10.0);
  PolyDG::Normal          n;
  PolyDG::Function        f(source), gd(uex);

  std::vector<PolyDG::BCLabelType> dirichlet = {1, 2, 3, 4, 5, 6};

  // Problem instantation and integration
  std::cout << "-------------------- Symmetric Problem --------------------" << std::endl;
  PolyDG::Problem poisson(Vh);

  poisson.integrateVol(dot(uGrad, vGrad), true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet, true);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  poisson.finalizeMatrix();

  Utilities::Watch ch;

  // LU
  std::cout << "\nSolving with SparseLU..." << std::endl;
  ch.start();
  poisson.solveLU();
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  ch.reset();

  // Chlolesky
  std::cout << "\nSolving with SparseCholesky..." << std::endl;
  ch.start();
  poisson.solveCholesky();
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  ch.reset();

  // Conjugate Gradient
  std::cout << "\nSolving with ConjugateGradient..." << std::endl;
  ch.start();
  poisson.solveCG(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim(), 1e-10);
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  ch.reset();

  // BiCGSTAB
  std::cout << "\nSolving with BiCGSTAB..." << std::endl;
  ch.start();
  poisson.solveBiCGSTAB(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim(), 1e-10);
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  // Another problem
  std::cout << "------------------ Non Symmetric Problem ------------------" << std::endl;
  poisson.clearMatrix();
  poisson.clearRhs();

  PolyDG::Stiff stiff;
  PolyDG::Mass  mass;

  poisson.integrateVol(stiff + mass, true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet,  false);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), false);

  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  poisson.finalizeMatrix();

  // LU
  std::cout << "\nSolving with SparseLU..." << std::endl;
  ch.start();
  poisson.solveLU();
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  ch.reset();

  // Chlolesky
  std::cout << "\nSolving with SparseCholesky..." << std::endl;
  ch.start();
  try
  {
    poisson.solveCholesky();
    ch.stop();
    std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
    std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
    std::cout << ch << std::endl;

  } catch(const std::exception&e)
  {
    std::cout << e.what() << std::endl;
  }

  ch.reset();

  // Conjugate Gradient
  std::cout << "\nSolving with ConjugateGradient..." << std::endl;
  ch.start();
  try
  {
    poisson.solveCG(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim(), 1e-10);
    ch.stop();
    std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
    std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
    std::cout << ch << std::endl;
  } catch(const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }

  ch.reset();

  // BiCGSTAB
  std::cout << "\nSolving with BiCGSTAB..." << std::endl;
  ch.start();
  poisson.solveBiCGSTAB(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim(), 1e-10);
  ch.stop();
  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;
  std::cout << ch << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;
}
