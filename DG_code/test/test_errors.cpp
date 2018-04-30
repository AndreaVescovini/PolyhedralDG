#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Problem.hpp"
#include "Operators.hpp"
#include "ExprOperators.hpp"
#include "Watch.hpp"

#include <Eigen/Core>

#include <cmath>

using namespace PolyDG;

int main()
{
  Timings::Watch ch;
  ch.start();

  std::vector<std::string> fileNames;
  // fileNames.push_back("../meshes/cube_str6t.mesh");
  // fileNames.push_back("../meshes/cube_str48t.mesh");
  // fileNames.push_back("../meshes/cube_str384t.mesh");
  // fileNames.push_back("../meshes/cube_str1296t.mesh");
  fileNames.push_back("/vagrant/pacs/progetto_codici/meshes/cube_str3072h.mesh");

  // auto uex = [](const Eigen::Vector3d& x) { return std::exp(x(0)*x(1)*x(2)); };
  // auto source = [&uex](const Eigen::Vector3d& x) { return -uex(x) * (x(0)*x(0)*x(1)*x(1) +
  //                                                                    x(1)*x(1)*x(2)*x(2) +
  //                                                                    x(0)*x(0)*x(2)*x(2) );};
  // auto uexGrad = [&uex](const Eigen::Vector3d& x) -> Eigen::Vector3d { return uex(x)*Eigen::Vector3d(x(1)*x(2),
  //                                                                                                    x(0)*x(2),
  //                                                                                                    x(0)*x(1)); };

  auto uex = [](Eigen::Vector3d x) -> double { return x(0); };
  auto source = [](Eigen::Vector3d /*x*/) -> double { return 0.0;};
  auto uexGrad = [](Eigen::Vector3d /*x*/)-> Eigen::Vector3d { return Eigen::Vector3d(1.0, 0.0, 0.0); };

  // auto uex = [](Eigen::Vector3d x) { return x(0)*x(0)*x(0)+10*x(1)*x(2)*x(2); };
  // auto source = [](Eigen::Vector3d x) { return -20*x(1)-6*x(0);};
  // auto uexGrad = [](Eigen::Vector3d x) { return Eigen::Vector3d(3*x(0)*x(0), 10*x(2)*x(2), 20*x(1)*x(2)); };

  PolyDG::MeshReaderPoly reader;

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

  unsigned r = 1;

  std::vector<double> errL2, errH10, hh;

  for(unsigned i = 0; i < fileNames.size(); i++)
  {
    PolyDG::Mesh Th(fileNames[i], reader);
    FeSpace Vh(Th, r);

    Problem prob(Vh);

    bool symform = true;
    prob.integrateVol(dot(uGrad, vGrad), symform);
    prob.integrateFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), symform);
    prob.integrateFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1, symform);
    prob.integrateVolRhs(f * v);
    prob.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, 1);

    prob.finalizeMatrix();

    ch.stop();
    prob.solveCG(Eigen::VectorXd::Zero(prob.getDim()), 1000);
    // prob.solveChol();
    // prob.solveLU();
    ch.start();

    errL2.push_back(prob.computeErrorL2(uex));
    errH10.push_back(prob.computeErrorH10(uexGrad));
    hh.push_back(Th.getMaxDiameter());

    std::cout << "Error L2 = " << errL2.back() << std::endl;
    std::cout << "Error H10 = " << errH10.back() << '\n' << std::endl;

    // std::cout << prob.getDim() << std::endl;
  }


  for(unsigned i = 0; i < fileNames.size() - 1; i++)
  {
    double orderConvL2 = std::log(errL2[i] / errL2[i+1]) / std::log(hh[i] / hh[i+1]);
    double orderConvH10 = std::log(errH10[i] / errH10[i+1]) / std::log(hh[i] / hh[i+1]);

    std::cout << "Order L2 = " << orderConvL2 << " - Order H10 = " << orderConvH10 << std::endl;
  }

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
