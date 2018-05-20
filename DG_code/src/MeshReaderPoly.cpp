/*!
    @file   MeshReaderPoly.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class MeshReaderPoly
*/

#include "MeshReaderPoly.hpp"
#include "FaceExt.hpp"
#include "PolyDG.hpp"
#include "Polyhedron.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <stdexcept>
#include <vector>

namespace PolyDG
{

MeshReaderPoly::MeshReaderPoly(const std::array<std::string, 4>& sections)
  : sections_{sections} {}

void MeshReaderPoly::read(Mesh& mesh, const std::string& fileName) const
{

  // I assume that entities are numerated from 1 to N in the meshFile,
  // then I store them in my mesh conteiners numerating from 0 to N-1.

  std::ifstream meshFile{fileName};
  if(meshFile.is_open() == false)
    throw std::runtime_error("Can't open mesh file. Mesh file does not exist or is corrupted.");

  // I use the proxy to access the Mesh class.
  MeshProxy mp(mesh);
  std::vector<Vertex>& vertList       = mp.getVerticesRef();
  std::vector<Tetrahedron>& tetraList = mp.getTetrahedraRef();
  std::vector<FaceExt>& faceExtList   = mp.getFacesExtRef();
  std::vector<Polyhedron>& polyList   = mp.getPolyhedraRef();

  bool found = goToSection(meshFile, 0);
  if(found == false)
    throw MeshFormatError("Error in mesh format, " + sections_[0] + " not found.");

  // Read vertices.
  SizeType verticesNo = 0;
  meshFile >> verticesNo;
  vertList.reserve(verticesNo);
  Vertex::resetCounter();

  std::array<Real, 3> curVertex;
  int label = 0;
  for(SizeType i = 0; i < verticesNo; i++)
  {
    meshFile >> curVertex[0] >> curVertex[1] >> curVertex[2] >> label;
    vertList.emplace_back(curVertex[0], curVertex[1], curVertex[2]);
  }

  found = goToSection(meshFile, 1);
  if(found == false)
    throw MeshFormatError("Error in mesh format, " + sections_[0] + " not found.");

  // Read tetrahedra.
  SizeType tetrahedraNo = 0;
  meshFile >> tetrahedraNo;
  tetraList.reserve(tetrahedraNo);
  Tetrahedron::resetCounter();

  std::array<unsigned, 4> curTet;
  for(SizeType i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> curTet[0] >> curTet[1] >> curTet[2] >> curTet[3] >> label;
    tetraList.emplace_back(vertList[curTet[0] - 1], vertList[curTet[1] - 1],
                           vertList[curTet[2] - 1], vertList[curTet[3] - 1]);
  }

  found = goToSection(meshFile, 2);
  if(found == false)
    throw MeshFormatError("Error in mesh format, " + sections_[0] + " not found.");

  // Read external faces.
  SizeType facesExtNo = 0;
  meshFile >> facesExtNo;
  faceExtList.reserve(facesExtNo);
  FaceExt::resetCounter();

  std::array<unsigned, 3> curFace;
  for(SizeType i = 0; i < facesExtNo; i++)
  {
    meshFile >> curFace[0] >> curFace[1] >> curFace[2] >> label;
    faceExtList.emplace_back(vertList[curFace[0] - 1], vertList[curFace[1] - 1],
                             vertList[curFace[2] - 1], label);
  }

  // If I don't find the section about polyhedra I consider every tetrahedron as
  // a a polyhedron.
  found = goToSection(meshFile, 3);

  // Read polyhedra.
  Polyhedron::resetCounter();
  SizeType polyhedraNo = 0;

  if(found == true)
  {
    meshFile >> polyhedraNo;
    polyList.resize(polyhedraNo);

    unsigned poly = 0;
    for(SizeType i = 0; i < tetrahedraNo; i++)
    {
      meshFile >> poly;
      polyList[poly - 1].addTetra(tetraList[i]);
      tetraList[i].setPoly(polyList[poly - 1]);
    }
  }
  else
  {
    polyhedraNo = tetrahedraNo;
    polyList.resize(polyhedraNo);

    for(SizeType i = 0; i < tetrahedraNo; i++)
    {
      polyList[i].addTetra(tetraList[i]);
      tetraList[i].setPoly(polyList[i]);
    }
  }

  meshFile.close();
}

bool MeshReaderPoly::goToSection(std::ifstream& meshFile, unsigned secNo) const
{
  std::string curLine = "";
  while(curLine != sections_[secNo] && curLine != "End")
    meshFile >> curLine;

  if(curLine == "End")
  {
    std::cout << sections_[3] << " not found." << std::endl;
    return false;
  }

  return true;
}

} // namespace PolyDG
