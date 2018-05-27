/*!
    @file   test_mesh.cpp
    @author Andrea Vescovini
    @brief  Test for the reading of the mesh
*/

#include "Mesh.hpp"
#include "MeshReaderPoly.hpp"
#include "Watch.hpp"

#include "GetPot.hpp"

#include <string>
#include <iostream>
#include <string>
#include <vector>

/*!
    All the meshes are read.
*/

int main(int argc, char* argv[])
{
  using PolyDG::Mesh;
  using PolyDG::MeshReaderPoly;

  Utilities::Watch ch;
  ch.start();

  GetPot comLine(argc, argv);
  const std::string fileName = comLine.follow("../data.pot", 2, "-f", "--file");
  GetPot fileData(fileName.c_str());

  const std::string meshDir = fileData("dir", "../../meshes");

  std::vector<std::string> fileNames;
  fileNames.emplace_back(meshDir + "/cube_str6t.mesh");
  fileNames.emplace_back(meshDir + "/cube_str48t.mesh");
  fileNames.emplace_back(meshDir + "/cube_str48p.mesh");
  fileNames.emplace_back(meshDir + "/cube_str48h.mesh");
  fileNames.emplace_back(meshDir + "/cube_str48ht.mesh");
  fileNames.emplace_back(meshDir + "/cube_str384t.mesh");
  fileNames.emplace_back(meshDir + "/cube_str384p.mesh");
  fileNames.emplace_back(meshDir + "/cube_str384h.mesh");
  fileNames.emplace_back(meshDir + "/cube_str384ht.mesh");
  fileNames.emplace_back(meshDir + "/cube_str1296t.mesh");
  fileNames.emplace_back(meshDir + "/cube_str1296p.mesh");
  fileNames.emplace_back(meshDir + "/cube_str1296h.mesh");
  fileNames.emplace_back(meshDir + "/cube_str1296ht.mesh");
  fileNames.emplace_back(meshDir + "/cube_str3072t.mesh");
  fileNames.emplace_back(meshDir + "/cube_str3072p.mesh");
  fileNames.emplace_back(meshDir + "/cube_str3072h.mesh");
  fileNames.emplace_back(meshDir + "/cube_str3072ht.mesh");

  MeshReaderPoly reader;

  Mesh Th(fileNames[0], reader);
  Th.printAll();

  Mesh Th2(fileNames[2], reader);
  Th2.printHead();
  Th2.exportMeshVTK("test_mesh.vtu");

  for(const std::string& f : fileNames)
  {
    Mesh Th3(f, reader);
    Th3.printInfo();
  }

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
