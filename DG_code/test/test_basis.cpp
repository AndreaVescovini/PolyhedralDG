// I test the initialization of the basis of the FeSpace.

#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"

#include <string>

int main()
{
  using PolyDG::FeSpace;
  using PolyDG::Mesh;
  using PolyDG::MeshReaderPoly;

  // Mesh Reading
  std::string fileName = "../meshes/cube_str48h.mesh";

  MeshReaderPoly reader;
  Mesh Th(fileName, reader);

  // FeSpace creation
  unsigned r = 1;
  FeSpace Vh(Th, r, 2, 2);

  // Print the computed values for the basis functions
  Vh.printInfo();
  Vh.printElemBasis();
  Vh.printElemBasisDer();
  Vh.printFaceBasis();
  Vh.printFaceBasisDer();

  return 0;
}
