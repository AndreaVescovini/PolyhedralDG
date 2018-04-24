#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Problem.hpp"
#include "Operators.hpp"
#include "ExprOperators.hpp"
#include "Watch.hpp"

#include <cmath>

using namespace PolyDG;

int main()
{
  std::string fileName = "/vagrant/pacs/progetto_codici/meshes/cube_str48p.mesh";

  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(fileName, reader);

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

  Function f([](Eigen::Vector3d /* x */) { return 0.0; });
  Function gd([](Eigen::Vector3d x) { return x(0); });

  Problem prob(Vh, true);

  bool symform = true;

  prob.integrateVol(dot(uGrad, vGrad), symform);
  prob.integrateFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), symform);
  prob.integrateFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1, symform);

  prob.integrateVolRhs(f * v);
  prob.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, 1);

  Timings::Watch ch;
  ch.start();

  std::cout << "\nSolving with SparseLU" << std::endl;
  prob.solveLU();
  // std::cout << prob.getSolution() << '\n' << std::endl;

  ch.stop();
  std::cout << ch << std::endl;
  ch.reset();
  ch.start();

  std::cout << "\nSolving with SparseCholesky" << std::endl;
  prob.solveChol();
  // std::cout << prob.getSolution() << '\n' <<  std::endl;

  ch.stop();
  std::cout << ch << std::endl;
  ch.reset();
  ch.start();

  std::cout << "\nSolving with ConjugateGradient" << std::endl;
  prob.solveCG(Eigen::VectorXd::Zero(prob.getDim()), 2 * prob.getDim());
  // std::cout << prob.getSolution() << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  // prob.exportSolutionVTK("output.vtu");

  return 0;
}
