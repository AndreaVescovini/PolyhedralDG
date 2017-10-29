#include "Mesh.hpp"

int main()
{
  geom::Mesh Th;
  Th.load("cube_str48h.mesh");
  Th.printHead();

  return 0;
}
