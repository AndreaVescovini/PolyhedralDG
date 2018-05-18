#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Watch.hpp"

#include <string>
#include <vector>

using PolyDG::Mesh;
using PolyDG::MeshReaderPoly;

int main()
{
  Timings::Watch ch;
  ch.start();

  std::vector<std::string> fileNames;
  fileNames.emplace_back("../meshes/cube_str6t.mesh");
  fileNames.emplace_back("../meshes/cube_str48t.mesh");
  fileNames.emplace_back("../meshes/cube_str48p.mesh");
  fileNames.emplace_back("../meshes/cube_str48h.mesh");
  fileNames.emplace_back("../meshes/cube_str48ht.mesh");
  fileNames.emplace_back("../meshes/cube_str384t.mesh");
  fileNames.emplace_back("../meshes/cube_str384p.mesh");
  fileNames.emplace_back("../meshes/cube_str384h.mesh");
  fileNames.emplace_back("../meshes/cube_str384ht.mesh");
  fileNames.emplace_back("../meshes/cube_str1296t.mesh");
  fileNames.emplace_back("../meshes/cube_str1296p.mesh");
  fileNames.emplace_back("../meshes/cube_str1296h.mesh");
  fileNames.emplace_back("../meshes/cube_str1296ht.mesh");
  fileNames.emplace_back("../meshes/cube_str3072t.mesh");
  fileNames.emplace_back("../meshes/cube_str3072p.mesh");
  fileNames.emplace_back("../meshes/cube_str3072h.mesh");
  fileNames.emplace_back("../meshes/cube_str3072ht.mesh");

  MeshReaderPoly reader;

  Mesh Th(fileNames[0], reader);
  Th.printAll();

  Mesh Th2(fileNames[2], reader);
  Th2.printHead();

  for(const std::string& f : fileNames)
  {
    Mesh Th3(f, reader);
    Th3.printInfo();
  }

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
