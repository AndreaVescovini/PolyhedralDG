#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Assembler.hpp"
#include "Operators.hpp"
#include "ExprOperators.hpp"

using namespace dgfem;

int main()
{
  std::string fileName = "../meshes/cube_str6t.mesh";

  geom::MeshReaderPoly reader;
  geom::Mesh Th(fileName, reader);

  unsigned r = 1;
  FeSpace Vh(Th, r);

  Stiff   stiff; // qui potrei assegnargli la viscosità forse, o forse no la viscosità dovrei semplicemente moltiplicargliela
  Mass    mass;
  GradPhiJ uGrad;
  GradPhiI vGrad;
  JumpPhiJ uJump;
  JumpPhiI vJump;
  AverGradPhiJ uGradAver;
  AverGradPhiI vGradAver;
  PenaltyScaling gamma(10.0);

  Assembler volumes(Vh);
  Assembler facesI(Vh);
  Assembler facesE(Vh);

  volumes.assembleVol(dot(uGrad, vGrad));
  facesI.assembleFacesInt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump));
  facesE.assembleFacesExt(-dot(uGradAver, vJump)-dot(uJump, vGradAver)+gamma*dot(uJump, vJump), 1);

  volumes.printMatrix();
  facesI.printMatrix();
  facesE.printMatrix();
  return 0;
}
