#include "Mesh.hpp"
#include "FeSpace.hpp"
#include "MeshReaderPoly.hpp"

int main()
{
  std::string fileName = "../meshes/cube_str6t.mesh";

  geom::MeshReaderPoly reader;
  geom::Mesh Th(fileName, reader);

  unsigned r = 1;
  dgfem::FeSpace Vh(Th, r, 2, 2);

  Vh.printElemBasis();
  Vh.printElemBasisDer();
  Vh.printFaceBasis();
  Vh.printFaceBasisDer();

  return 0;
}
