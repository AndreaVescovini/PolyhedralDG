// I test the initialization of the basis of the FeSpace.

#include "FeSpace.hpp"
#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Watch.hpp"

#include "GetPot.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
  using PolyDG::FeSpace;
  using PolyDG::Mesh;
  using PolyDG::MeshReaderPoly;

  Utilities::Watch ch;
  ch.start();

  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshFile = fileData("dir", "../../meshes") + "/cube_str48h.mesh";

  // Mesh Reading
  MeshReaderPoly reader;
  Mesh Th(meshFile, reader);

  // FeSpace creation
  unsigned r = 1;
  FeSpace Vh(Th, r, 2, 2);

  // Print the computed values for the basis functions
  Vh.printInfo();
  Vh.printElemBasis();
  Vh.printElemBasisDer();
  Vh.printFaceBasis();
  Vh.printFaceBasisDer();

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
