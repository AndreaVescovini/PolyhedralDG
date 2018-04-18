#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "FeSpace.hpp"

#include <string>

int main()
{
  using PolyDG::MeshReaderPoly;
  using PolyDG::Mesh;
  using PolyDG::FeSpace;

  std::string fileName = "../meshes/cube_str6t.mesh";

  MeshReaderPoly reader;
  Mesh Th(fileName, reader);
  Th.printHead();
  // Th.printAll();

  unsigned r = 1;
  FeSpace Vh(Th, r);


  return 0;
}
