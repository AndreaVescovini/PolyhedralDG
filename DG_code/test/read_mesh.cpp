#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"

int main()
{
  using dgfem::MeshReaderPoly;
  using dgfem::Mesh;

  MeshReaderPoly reader({"Vertices", "Tetrahedra", "Polyhedra"});
  Mesh Th("../meshes/cube_str48h.mesh", reader);
  Th.printHead();

  return 0;
}
