#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"
#include "Assembler.hpp"
#include "Expr.hpp"
#include "Operators.hpp"

using namespace dgfem;

int main()
{
  std::string fileName = "../meshes/cube_str6t.mesh";

  geom::MeshReaderPoly reader;
  geom::Mesh Th(fileName, reader);

  unsigned r = 2;
  FeSpace Vh(Th, r);

  Expr<Stiff> stiff;
  Assembler<Stiff> exprTemplate(stiff, Vh);
  exprTemplate.assembleVol();
  exprTemplate.printMatrix();

  return 0;
}
