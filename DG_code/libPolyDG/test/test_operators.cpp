/*!
    @file   test_operators.cpp
    @author Andrea Vescovini
    @brief  Test for different operators
*/

#include "ExprOperators.hpp"
#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Problem.hpp"
#include "Utilities.hpp"
#include "Watch.hpp"

#include <Eigen/Core>
#include "GetPot.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

/*!
    Three different problems are tested with a konwn polynomial solution:

    1) \f$ - \Delta u = f_1  \quad \text{in} \quad \Omega\\
                    u = g_d  \quad \text{on} \quad \partial \Omega \f$

    2)  \f$ - \Delta u + div(\mathbf{b}u) = f_2  \quad \text{in} \quad \Omega\\
                                        u = g_d  \quad \text{on} \quad \partial \Omega \f$

    3)           \f$ - \Delta u = f_3 \quad \text{in} \quad \Omega\\
                              u = g_d \quad \text{on} \quad \Gamma_{Dirichlet}\\
      \nabla u \cdot \mathbf{n} = g_n \quad \text{on} \quad \Gamma_{Neumann} \f$

*/

int main(int argc, char* argv[])
{
  using Utilities::pow;

  Utilities::Watch ch;
  ch.start();

  // Exact solution
  auto uex     = [](const Eigen::Vector3d& x) { return x(0) * x(1); };
  auto uexGrad = [](const Eigen::Vector3d& x) { return Eigen::Vector3d(x(1), x(0), 0.0); };

  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshFile = fileData("dir", "../../meshes") + "/cube_str384ht.mesh";

  // Mesh reading
  PolyDG::MeshReaderPoly reader;
  PolyDG::Mesh Th(meshFile, reader);
  Th.printInfo();

  // Operators
  PolyDG::Stiff           stiff;
  PolyDG::Mass            mass;
  PolyDG::PhiJ            u;
  PolyDG::PhiI            v;
  PolyDG::GradPhiJ        uGrad;
  PolyDG::GradPhiI        vGrad;
  PolyDG::JumpPhiJ        uJump;
  PolyDG::JumpPhiI        vJump;
  PolyDG::AverPhiJ        uAver;
  PolyDG::AverGradPhiJ    uGradAver;
  PolyDG::AverGradPhiI    vGradAver;
  PolyDG::PenaltyScaling  gamma(10.0);
  PolyDG::Normal          n;
  PolyDG::Function        gd(uex);

  std::cout << "---------------- Diffusion-reaction Problem ---------------" << std::endl;

  // Source term
  PolyDG::Function f1([](const Eigen::Vector3d& x) { return x(0) * x(1); });

  // FeSpace Creation
  PolyDG::FeSpace Vh1(Th, 2);

  // Problem instantation and integration
  PolyDG::Problem diffreac(Vh1);

  std::vector<PolyDG::BCLabelType> dirichlet = {1, 2, 3, 4, 5, 6};

  diffreac.integrateVol(stiff + mass, true);
  diffreac.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet, true);
  diffreac.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
  diffreac.integrateVolRhs(f1 * v);
  diffreac.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  diffreac.finalizeMatrix();

  // Sparse Chlolesky
  std::cout << "Solving with SparseCholesky..." << std::endl;
  diffreac.solveCholesky();

  std::cout << "L2  error = " << diffreac.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << diffreac.computeErrorH10(uexGrad) << std::endl;

  std::cout << "----------- Advection-diffusion-reaction Problem ----------" << std::endl;

  // Source term and advective field
  PolyDG::Function3 b([](const Eigen::Vector3d& x) { return Eigen::Vector3d(x(0), 0.0, x(2)); });
  PolyDG::Function f2([](const Eigen::Vector3d& x) { return 4.0 * x(0) * x(1); });

  // FeSpace Creation
  unsigned r = 2;
  PolyDG::FeSpace Vh2(Th, r, 2 * r, 2 * r + 1);

  // Problem instantation and integration
  PolyDG::Problem adr(Vh2);

  adr.integrateVol(stiff + mass - dot(b * u, vGrad), false);
  adr.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump) + dot(b * u, n) * v, dirichlet, true);
  adr.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump) + dot(b, vJump) * uAver, false);
  adr.integrateVolRhs(f2 * v);
  adr.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  adr.finalizeMatrix();

  // Sparse LU
  std::cout << "Solving with SparseLU..." << std::endl;
  adr.solveLU();

  std::cout << "L2  error = " << adr.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << adr.computeErrorH10(uexGrad) << std::endl;

  std::cout << "------ Poisson problem with mixed boundary conditions -----" << std::endl;

  // Source term and Neumann b.c.
  PolyDG::Function f3([](const Eigen::Vector3d& /* x */) { return 0.0; });
  PolyDG::Function  h([](const Eigen::Vector3d& x) { return x(1) * (2.0 * x(0) - 1.0); });

  // Problem instantation and integration
  PolyDG::Problem poisson(Vh1);

  dirichlet = {3, 4, 5, 6};
  std::vector<PolyDG::BCLabelType> neumann = {1, 2};

  poisson.integrateVol(stiff, true);
  poisson.integrateFacesExt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), dirichlet, true);
  poisson.integrateFacesInt(-dot(uGradAver, vJump) - dot(uJump, vGradAver) + gamma * dot(uJump, vJump), true);
  poisson.integrateVolRhs(f3 * v);
  poisson.integrateFacesExtRhs(h * v, neumann);
  poisson.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v, dirichlet);

  poisson.finalizeMatrix();

  // Sparse LU
  std::cout << "Solving with SparseCholesky..." << std::endl;
  poisson.solveCholesky();

  std::cout << "L2  error = " << poisson.computeErrorL2(uex) << std::endl;
  std::cout << "H10 error = " << poisson.computeErrorH10(uexGrad) << std::endl;

  std::cout << "-----------------------------------------------------------" << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
