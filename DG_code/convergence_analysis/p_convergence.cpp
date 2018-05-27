// This example computes the p-convergence order of the method solving the Poisson
// problem:
//      -laplacian(u) = f   in omega
//                 u  = gd  in delta_omega
//
// employing each time a FeSpace of higher degree. The data f and gd are
// chosen such that the exact solution u is u = exp(x * y * z).

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
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
  using namespace PolyDG;
  using Utilities::pow;

  Utilities::Watch ch;
  ch.start();

  // Mesh file
  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("data_pconv.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshFile = fileData("directory", "../../meshes") + '/' + fileData("mesh", "cube_str1296p.mesh");

  // Maximum degree of the FeSpace
  const unsigned degMax = fileData("degMax", 5);

  // Solver
  const std::string solverType = fileData("solverType", "iterative");
  const unsigned maxIter       = fileData("maxIter", 10000);
  const double tolerance       = fileData("tolerance", 1e-10);

  // Export the solutions
  const std::string exportSol = fileData("exportSol", "no");

  // Functions
  auto uex = [](const Eigen::Vector3d& x) { return std::exp(x(0) * x(1) * x(2)); };
  auto source = [&uex](const Eigen::Vector3d& x) {
    return -uex(x) * (pow(x(0) * x(1), 2) + pow(x(1) * x(2), 2) + pow(x(0) * x(2), 2));};
  auto uexGrad = [&uex](const Eigen::Vector3d& x) -> Eigen::Vector3d {
    return uex(x) * Eigen::Vector3d(x(1) * x(2),  x(0) * x(2), x(0) * x(1)); };

  // Read the mesh
  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(meshFile, reader);

  // Operators
  PhiI v;
  GradPhiJ uGrad;
  GradPhiI vGrad;
  JumpPhiJ uJump;
  JumpPhiI vJump;
  AverGradPhiJ uGradAver;
  AverGradPhiI vGradAver;
  PenaltyScaling gamma(10.0);
  Normal n;
  Function f(source);
  Function gd(uex);

  // Boundary labels
  std::vector<PolyDG::BCLabelType> dirichlet = {1, 2, 3, 4, 5, 6};

  std::vector<double> errL2, errH10, orderConvL2, orderConvH10;

  for(unsigned deg = 1; deg <= degMax; deg++)
  {
    // FeSpace creation
    FeSpace Vh(Th, deg);
    Vh.printInfo();

    // Problem instantation
    Problem poisson(Vh);

    // Integration of the variational form (symmetric)
    poisson.integrateVol(dot(uGrad, vGrad), true);
    poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
    poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet, true);
    poisson.integrateVolRhs(f * v);
    poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

    poisson.finalizeMatrix();

    // Solve the linear system
    if(solverType == "direct")
      poisson.solveCholesky();
    else
      poisson.solveCG(Eigen::VectorXd::Zero(poisson.getDim()), maxIter, tolerance);

    // Save the solution
    if(exportSol == "yes")
    {
      std::string solName = "solution_" + std::to_string(deg) + ".vtu";
      poisson.exportSolutionVTK(solName);
    }

    // Computation of the errors
    errL2.push_back(poisson.computeErrorL2(uex));
    errH10.push_back(poisson.computeErrorH10(uexGrad));

    std::cout << "Error L2  = " << errL2.back() << std::endl;
    std::cout << "Error H10 = " << errH10.back() << '\n' << std::endl;

    if(deg >= 2)
    {
      orderConvL2.push_back(std::log(errL2[deg - 2] / errL2[deg - 1]) /
                            std::log(static_cast<double>(deg) / (deg - 1)));
      orderConvH10.push_back(std::log(errH10[deg - 2] / errH10[deg - 1]) /
                             std::log(static_cast<double>(deg) / (deg - 1)));
    }
  }

  // Print the results
  std::cout << "Order L2  =";
  for(const double& i : orderConvL2)
    std::cout << ' ' << i;

  std::cout << "\nOrder H10 =";
  for(const double& i : orderConvH10)
    std::cout << ' ' << i;

  std::cout << '\n' << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
