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
  FeSpace Vh(Th, r, 2, 2);

  Stiff   stiff; // qui potrei assegnargli la viscosità forse, o forse no la viscosità dovrei semplicemente moltiplicargliela
  Mass    mass;
  PhiI v;
  GradPhiJ uGrad;
  GradPhiI vGrad;
  JumpPhiJ uJump;
  JumpPhiI vJump;
  AverGradPhiJ uGradAver;
  AverGradPhiI vGradAver;
  PenaltyScaling gamma(10.0);
  Normal n;

  auto source = [](Eigen::Vector3d x) {
    return -std::exp(x(0)*x(1)*x(2)) * (x(0)*x(0)*x(1)*x(1) +
                                        x(1)*x(1)*x(2)*x(2) +
                                        x(0)*x(0)*x(2)*x(2) ); };

  auto dirichlet = [](Eigen::Vector3d x) { return std::exp(x(0)*x(1)*x(2)); };

  Function f(source);
  Function gd(dirichlet);

  Assembler volumes(Vh);
  Assembler facesI(Vh);
  Assembler facesE(Vh);

  // volumes.assembleVol(dot(uGrad, vGrad));
  // facesI.assembleFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump));
  // facesE.assembleFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1);

  // volumes.assembleVolRhs(source * v);
  facesE.assembleFacesExtRhs(-gd * dot(n, vGrad) + gamma * gd * v);

  // volumes.printRhs();
  facesE.printRhs();
  // volumes.printMatrix();
  // facesI.printMatrix();
  // facesE.printMatrix();
  return 0;
}
