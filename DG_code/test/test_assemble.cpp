#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Assembler.hpp"
#include "Operators.hpp"

using namespace dgfem;

int main()
{
  std::string fileName = "../meshes/cube_str6t.mesh";

  geom::MeshReaderPoly reader;
  geom::Mesh Th(fileName, reader);

  unsigned r = 1;
  FeSpace Vh(Th, r);

  // Expr<Stiff> stiff;
  // Assembler<Stiff> exprTemplate(stiff, Vh);
  // exprTemplate.assembleVol();

  Stiff stiff; // qui potrei assgnarli la viscosità forse, o forse no la viscosità dovrei semplicemente moltiplicargliela
  Assembler exprTemplate(Vh);
  exprTemplate.assembleVol(stiff);

  exprTemplate.printMatrix();

  return 0;
}
