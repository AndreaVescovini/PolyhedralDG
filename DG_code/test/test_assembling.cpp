#include "ExprOperators.hpp"
#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Operators.hpp"
#include "Problem.hpp"

#include <cmath>

int main()
{
  // Mesh reading
  std::string fileName = "../meshes/cube_str48h.mesh";

  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(fileName, reader);

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
    return -std::exp(x(0) * x(1) * x(2)) * (x(0) * x(0) * x(1) * x(1) +
                                            x(1) * x(1) * x(2) * x(2) +
                                            x(0) * x(0) * x(2) * x(2)); });
  PolyDG::Function gd([](Eigen::Vector3d x) { return std::exp(x(0) * x(1) * x(2)); });

  // Problem instantation and integration
  PolyDG::Problem poisson(Vh);

  poisson.integrateVol(dot(uGrad, vGrad), true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), 1,  true);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);

  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, 1);

  poisson.finalizeMatrix();

  poisson.printInfo();
  std::cout << poisson.getMatrix() << std::endl;
  std::cout << poisson.getRhs() << std::endl;

  // Another poroblem
  poisson.clearMatrix();
  poisson.clearRhs();

  PolyDG::Stiff stiff;
  PolyDG::Mass  mass;

  poisson.integrateVol(stiff + mass, true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), 1,  false);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) + dot(uJump, vGradAver) + gamma * dot(uJump, vJump), false);

  poisson.integrateVolRhs(f * v);
  poisson.integrateFacesExtRhs(gd * dot(n, vGrad) + gamma * gd * v, 1);

  poisson.finalizeMatrix();

  poisson.printInfo();
  std::cout << poisson.getMatrix() << std::endl;
  std::cout << poisson.getRhs() << std::endl;

  return 0;
}
