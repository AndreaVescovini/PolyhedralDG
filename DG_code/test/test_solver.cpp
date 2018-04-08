#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Assembler.hpp"
#include "Operators.hpp"
#include "ExprOperators.hpp"
#include <cmath>

using namespace dgfem;

int main()
{
  std::string fileName = "../meshes/cube_str6t.mesh";

  geom::MeshReaderPoly reader;
  geom::Mesh Th(fileName, reader);

  unsigned r = 1;
  FeSpace Vh(Th, r);

  PhiI v;
  GradPhiJ uGrad;
  GradPhiI vGrad;
  JumpPhiJ uJump;
  JumpPhiI vJump;
  AverGradPhiJ uGradAver;
  AverGradPhiI vGradAver;
  PenaltyScaling gamma(10.0);
  Normal n;

  Function f([](Eigen::Vector3d x) { return 0.0; });
  Function gd([](Eigen::Vector3d x) { return x(0); });

  Assembler prob(Vh, true);

  bool symform = true;

  prob.integrateVol(dot(uGrad, vGrad), symform);
  prob.integrateFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), symform);
  prob.integrateFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1, symform);

  prob.integrateVolRhs(f * v);
  prob.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v);

  std::cout << "\nSolving with SparseLU" << std::endl;
  prob.solveLU();
  std::cout << prob.getSolution() << '\n' << std::endl;

  std::cout << "\nSolving with SparseCholesky" << std::endl;
  prob.solveChol();
  std::cout << prob.getSolution() << '\n' <<  std::endl;

  std::cout << "\nSolving with ConjugateGradient" << std::endl;
  prob.solveCG(Eigen::VectorXd::Zero(prob.getDim()), 100);
  std::cout << prob.getSolution() << std::endl;

  return 0;
}
