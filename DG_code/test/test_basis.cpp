#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"

#include <string>

using PolyDG::FeSpace;
using PolyDG::Mesh;
using PolyDG::MeshReaderPoly;

int main()
{
  std::string fileName = "../meshes/cube_str48h.mesh";

  MeshReaderPoly reader;
  Mesh Th(fileName, reader);

  unsigned r = 1;
  FeSpace Vh(Th, r, 2, 2);

  Vh.printInfo();

  Vh.printElemBasis();
  Vh.printElemBasisDer();
  Vh.printFaceBasis();
  Vh.printFaceBasisDer();

  return 0;
}
