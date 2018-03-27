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

  // for(int i = 0; i < Th.getFacesIntNo(); i++)
  // {
  //   auto faccia = Th.getFaceInt(i);
  //   std::cout << "Normale" << '\n';
  //   std::cout << faccia.getNormal().transpose() << std::endl;
  //   std::cout << "Vertici faccia" << '\n';
  //   for(int j = 0; j < 3; j++)
  //     std::cout << faccia.getVertex(j).getCoords().transpose() << std::endl;
  //   std::cout << "Vertici Tetra 1 e poi 2" << '\n';
  //   std::cout << faccia.getTet1().getVertex(3-faccia.getFaceNoTet1()).getCoords().transpose() << std::endl;
  //   std::cout << faccia.getFaceNoTet2() << std::endl;
  //
  //
  //   std::cout << "--------------------------------------------" << '\n';
  // }

  unsigned r = 1;
  FeSpace Vh(Th, r);

  Stiff   stiff; // qui potrei assegnargli la viscosità forse, o forse no la viscosità dovrei semplicemente moltiplicargliela
  Mass    mass;
  GradPhi gradPhi;
  PhiI    u;
  PhiJ    v;
  Phi     phi;
  JumpIntPhi jumpIntPhi;
  AverIntGradPhi averIntGradPhi;
  PenaltyScaling gamma(10.0);

  Assembler volumes(Vh);
  Assembler facesI(Vh);
  Assembler facesIt(Vh);
  // volumes.assembleVol(dot(gradPhi, gradPhi));

  facesI.assembleFacesInt(-dot(jumpIntPhi, averIntGradPhi)-dot(averIntGradPhi, jumpIntPhi), true);
  // facesIt.assembleFacesInt();

  facesI.printMatrix();
  // facesIt.printMatrix();
  return 0;
}
