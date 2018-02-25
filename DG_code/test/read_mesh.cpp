#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"

int main()
{
  using dgfem::MeshReaderPoly;
  using dgfem::Mesh;

  MeshReaderPoly reader;
  Mesh Th("../meshes/cube_str6t.mesh", reader);
  Th.printHead();
  // Th.printAll();

  return 0;
}
