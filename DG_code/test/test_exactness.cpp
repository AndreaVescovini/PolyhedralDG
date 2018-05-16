#include "ExprOperators.hpp"
#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Operators.hpp"
#include "Problem.hpp"

#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <vector>

int main()
{
  // Mesh reading
  std::string fileName = "../meshes/cube_str384h.mesh";

  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(fileName, reader);
  Th.printInfo();

  // Function definition
  std::vector<std::function<PolyDG::Real (const Eigen::Vector3d&)>> uex;
  std::vector<std::function<PolyDG::Real (const Eigen::Vector3d&)>> source;
  std::vector<std::function<Eigen::Vector3d (const Eigen::Vector3d&)>> uexGrad;

  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(1); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(0) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(0) * x(0) * x(1) * x(2); });
  uex.emplace_back([](const Eigen::Vector3d& x) { return x(0) * x(0) * x(0) * x(1) * x(1) * x(2); });

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
                                                                             x(0) * x(0) * x(0) * x(2),
                                                                             x(0) * x(0) * x(0) * x(1)); });
  uexGrad.emplace_back([](const Eigen::Vector3d& x) { return Eigen::Vector3d(3.0 * x(0) * x(0) * x(1) * x(1) * x(2),
                                                                             2.0 * x(0) * x(0) * x(0) * x(1) * x(2),
                                                                             x(0) * x(0) * x(0) * x(1) * x(1)); });
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

  std::vector<double> errL2, errH10, hh;

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
    poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), 1, true);
    poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
    poisson.integrateVolRhs(f * v);
    poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, 1);

    poisson.finalizeMatrix();

    poisson.solveCG(Eigen::VectorXd::Zero(poisson.getDim()), 2 * poisson.getDim());

    std::cout << "Error L2 =  " << poisson.computeErrorL2(uex[r - 1]) << std::endl;
    std::cout << "Error H10 = " << poisson.computeErrorH10(uexGrad[r - 1]) << '\n' << std::endl;
  }

  return 0;
}