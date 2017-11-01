#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"

int main()
{
  dgfem::MeshReaderPoly reader("../meshes/", "Vertices", "Tetrahedra", "Polyhedra");
  dgfem::Mesh Th("cube_str48h.mesh", reader);
  Th.printHead();

  return 0;
}
