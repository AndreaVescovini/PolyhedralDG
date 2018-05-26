// Test for the assembling of the matrix of the linear system.
//
// 1)  -laplacian(u) = f   in omega
//                u = gd   in delta_omega
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

#include "GetPot.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
  using Utilities::pow;

  Utilities::Watch ch;
  ch.start();

  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshFile = fileData("dir", "../../meshes") + "/cube_str48h.mesh";

  // Mesh reading
  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(meshFile, reader);

  // FeSpace creation
  unsigned r = 1;
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

  PolyDG::Function f([](Eigen::Vector3d x) {
    return -std::exp(x(0) * x(1) * x(2)) * (pow(x(0) * x(1), 2) + pow(x(1) * x(2), 2) + pow(x(0) * x(2), 2)); });
  PolyDG::Function gd([](Eigen::Vector3d x) { return std::exp(x(0) * x(1) * x(2)); });

  // Problem instantation and integration
  PolyDG::Problem poisson(Vh);

  std::vector<PolyDG::BCLabelType> dirichlet = {1, 2, 3, 4, 5, 6};

  poisson.integrateVol(dot(uGrad, vGrad), true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet,  true);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);

  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  poisson.finalizeMatrix();

  poisson.printInfo();
  std::cout << "Matrix:\n" <<  poisson.getMatrix() << std::endl;
  std::cout << "Rhs:\n" << poisson.getRhs() << std::endl;

  // Another poroblem.
  poisson.clearMatrix();
  poisson.clearRhs();

  PolyDG::Stiff stiff;
  PolyDG::Mass  mass;
  PolyDG::PhiJ  u;

  poisson.integrateVol(stiff + u * v, true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet,  false);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), false);

  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  poisson.finalizeMatrix();

  poisson.printInfo();
  std::cout << "Matrix:\n" <<  poisson.getMatrix() << std::endl;
  std::cout << "Rhs:\n" << poisson.getRhs() << '\n' << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
