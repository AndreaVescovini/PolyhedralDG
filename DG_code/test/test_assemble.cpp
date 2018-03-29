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
  GradPhi gradPhi;
  PhiI    u;
  PhiJ    v;
  Phi     phi;
  JumpPhi jumpIntPhi;
  AverGradPhi averIntGradPhi;
  PenaltyScaling gamma(10.0);

  Assembler volumes(Vh);
  Assembler facesI(Vh);
  Assembler facesE(Vh);

  // volumes.assembleVol(dot(gradPhi, gradPhi));
  // facesI.assembleFacesInt(-dot(averIntGradPhi, jumpIntPhi)-dot(jumpIntPhi, averIntGradPhi)+gamma*dot(jumpIntPhi,jumpIntPhi));
  facesE.assembleFacesExt(-dot(averIntGradPhi, jumpIntPhi)-dot(jumpIntPhi, averIntGradPhi)+gamma*dot(jumpIntPhi,jumpIntPhi), 1);

  facesE.printMatrix();
  // facesIt.printMatrix();
  return 0;
}
