/*!
    @file   test_exactness.cpp
    @author Andrea Vescovini
    @brief  Test for the exactness of the method
*/

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
#include <functional>
#include <iostream>
#include <string>
#include <vector>

/*!
    The exactness of the method applied to problems with polynomial solutions is tested,
    the error is expected to be 0 if the order of the FeSpace is greater of equal to
    that of the solution and the integration is exact.@n
    The problem is:

    \f$ - \Delta u = f  \quad \text{in} \quad \Omega\\
                     u = g_d  \quad \text{on} \quad \partial \Omega \f$
*/

int main(int argc, char* argv[])
{
  Utilities::Watch ch;
  ch.start();

  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshFile = fileData("dir", "../../meshes") + "/cube_str384h.mesh";

  // Mesh reading
  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(meshFile, reader);
  Th.printInfo();

  // Function definition
  std::vector<std::function<PolyDG::Real (const Eigen::Vector3d&)>> uex;
  std::vector<std::function<PolyDG::Real (const Eigen::Vector3d&)>> source;
  std::vector<std::function<Eigen::Vector3d (const Eigen::Vector3d&)>> uexGrad;

  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(1); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(0) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return pow(x(0), 3) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return pow(x(0), 3) * x(1) * x(1) * x(2); });

  source.emplace_back([](const Eigen::Vector3d& /* x */) { return 0.0; });
  source.emplace_back([](const Eigen::Vector3d& /* x */) { return 0.0; });
  source.emplace_back([](const Eigen::Vector3d& /* x */) { return 0.0; });
  source.emplace_back([](const Eigen::Vector3d& x) { return -2.0 * x(1) * x(2); });
  source.emplace_back([](const Eigen::Vector3d& x) { return -6.0 * x(0) * x(1) * x(2); });
  source.emplace_back([](const Eigen::Vector3d& x) { return -2.0 * x(0) * x(2) * (3.0 * x(1) * x(1) + x(0) * x(0)); });

  uexGrad.emplace_back([](const Eigen::Vector3d& /* x */) { return Eigen::Vector3d(1.0,0.0, 0.0); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(x(1), x(0), 0.0); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(x(1) * x(2),
                                                                             x(0) * x(2),
                                                                             x(0) * x(1)); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(2.0 * x(0) * x(1) * x(2),
                                                                             x(0) * x(0) * x(2),
                                                                             x(0) * x(0) * x(1)); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(3.0 * x(0) * x(0) * x(1) * x(2),
                                                                             pow(x(0), 3) * x(2),
                                                                             pow(x(0), 3) * x(1)); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(3.0 * pow(x(0) * x(1), 2) * x(2),
                                                                             2.0 * pow(x(0), 3) * x(1) * x(2),
                                                                             x(0) * pow(x(0) * x(1), 2)); });
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

  std::vector<PolyDG::BCLabelType> dirichlet = {1, 2, 3, 4, 5, 6};

  std::vector<double> errL2, errH10, hh;

  std::cout << "Testing the exactness of the method..." << std::endl;

  for(unsigned r = 1; r < 7; r++)
  {
    std::cout << "Degree = " << r << std::endl;

    PolyDG::Function f(source[r - 1]);
    PolyDG::Function gd(uex[r - 1]);

    // FeSpace creation
    PolyDG::FeSpace Vh(Th, r);

    // Problem instantation and integration
    PolyDG::Problem poisson(Vh);

    poisson.integrateVol(dot(uGrad, vGrad), true);
    poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet, true);
    poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
    poisson.integrateVolRhs(f * v);
    poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

    poisson.finalizeMatrix();

    poisson.solveCG(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim());

    std::cout << "Error L2  = " << poisson.computeErrorL2(uex[r - 1]) << std::endl;
    std::cout << "Error H10 = " << poisson.computeErrorH10(uexGrad[r - 1]) << '\n' << std::endl;
  }

  std::cout << "Test finished" << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
