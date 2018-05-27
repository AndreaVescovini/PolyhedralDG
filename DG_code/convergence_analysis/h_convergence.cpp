// This example computes the h-convergence order of the method solving the Poisson
// problem:
//      -laplacian(u) = f   in omega
//                 u  = gd  in delta_omega
//
// over a sequence of four meshes more and more refined. The data f and gd are
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
  const std::string fileName = comLine.follow("data_hconv.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshDir = fileData("directory", "../../meshes");

  std::vector<std::string> fileNames;
  fileNames.push_back(meshDir + "/cube_str48"   + fileData("meshType", "h") + ".mesh");
  fileNames.push_back(meshDir + "/cube_str384"  + fileData("meshType", "h") + ".mesh");
  fileNames.push_back(meshDir + "/cube_str1296" + fileData("meshType", "h") + ".mesh");
  fileNames.push_back(meshDir + "/cube_str3072" + fileData("meshType", "h") + ".mesh");

  // Degree of the FeSapce
  const unsigned deg = fileData("degree", 2);

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

  PolyDG::MeshReaderPoly reader;

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

  std::vector<double> errL2, errH10, hh, orderConvL2, orderConvH10;

  for(SizeType i = 0; i < fileNames.size(); i++)
  {
    // Reading the mesh
    PolyDG::Mesh Th(fileNames[i], reader);
    Th.printInfo();

    // FeSpace creation
    FeSpace Vh(Th, deg);

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
      std::string solName = "solution_" + std::to_string(Th.getTetrahedraNo()) + ".vtu";
      poisson.exportSolutionVTK(solName);
    }

    // Computation of the errors
    errL2.push_back(poisson.computeErrorL2(uex));
    errH10.push_back(poisson.computeErrorH10(uexGrad));
    hh.push_back(Th.getMaxDiameter());

    std::cout << "Error L2  = " << errL2.back() << std::endl;
    std::cout << "Error H10 = " << errH10.back() << '\n' << std::endl;

    if(i > 0)
    {
      orderConvL2.push_back(std::log(errL2[i - 1] / errL2[i]) / std::log(hh[i - 1] / hh[i]));
      orderConvH10.push_back(std::log(errH10[i - 1] / errH10[i]) / std::log(hh[i - 1] / hh[i]));
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
