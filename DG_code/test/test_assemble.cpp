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

  Function f([](Eigen::Vector3d x) {
    return -std::exp(x(0)*x(1)*x(2)) * (x(0)*x(0)*x(1)*x(1) +
                                        x(1)*x(1)*x(2)*x(2) +
                                        x(0)*x(0)*x(2)*x(2) ); });
  Function gd([](Eigen::Vector3d x) { return std::exp(x(0)*x(1)*x(2)); });

  Assembler volumes(Vh);
  Assembler facesI(Vh);
  Assembler facesE(Vh);

  // volumes.integrateVol(dot(uGrad, vGrad), true);
  facesI.integrateFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), true);
  // facesE.integrateFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1, true);

  // volumes.integrateVolRhs(f * v);
  // volumes.integrateFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v);

  // volumes.printRhs();
  // volumes.printMatrix();
  facesI.printMatrixSym();
  // facesE.printMatrix();
  return 0;
}
