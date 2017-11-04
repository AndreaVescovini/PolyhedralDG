#include "MeshReaderPoly.hpp"

#include <vector>
#include <fstream>
#include "Point.hpp"
#include "Tetrahedron.hpp"
#include "FaceExt.hpp"

namespace dgfem {

MeshReaderPoly::MeshReaderPoly(std::array<std::string, 4> sections)
  : sections_{sections} {}

void MeshReaderPoly::read(Mesh& mesh, const std::string& fileName) const
{
  using geom::Point;
  using geom::Tetrahedron;
  using geom::FaceExt;
  using geom::real;
  using geom::labelType;

  std::ifstream meshFile{fileName};
  if( meshFile.is_open() == false )
    std::cerr << "Can't open mesh file. Mesh file does not exist or is corrupted" << std::endl;

  // Use the proxy to access the Mesh class
  MeshProxy mp(mesh);
  std::vector<Point>& vertList = mp.getVerticesRef();
  std::vector<Tetrahedron>& tetraList = mp.getTetrahedraRef();
  std::vector<FaceExt>& faceExtList = mp.getFacesExtRef();

  goToSection(meshFile, 0);

  // Read vertices
  size_t verticesNo = 0;
  meshFile >> verticesNo;
  vertList.reserve(verticesNo);

  std::array<real, 3> curVertex;
  unsigned label = 0; // used to read unuseful labels
  for(size_t i = 0; i < verticesNo; i++)
  {
    meshFile >> curVertex[0] >> curVertex[1] >> curVertex[2] >> label;
    vertList.emplace_back(curVertex);
  }

  goToSection(meshFile, 1);

  // Read tetrahedra
  size_t tetrahedraNo = 0;
  meshFile >> tetrahedraNo;
  tetraList.reserve(tetrahedraNo);

  std::array<labelType, 4> curTetra;
  for(size_t i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> curTetra[0] >> curTetra[1] >> curTetra[2] >> curTetra[3] >> label;
    tetraList.emplace_back(curTetra);
  }

  goToSection(meshFile, 2);

  // Read Faces
  size_t facesExtNo = 0;
  meshFile >> facesExtNo;
  faceExtList.reserve(facesExtNo);

  std::array<labelType, 3> curFace;
  for(size_t i = 0; i < facesExtNo; i++)
  {
    meshFile >> curFace[0] >> curFace[1] >> curFace[2] >> label;
    faceExtList.emplace_back(curFace, label);
  }

  goToSection(meshFile, 3);

  // Read polyhedra
  size_t polyhedraNo = 0;
  meshFile >> polyhedraNo;
  // polyList.reserve(polyhedraNo);

  labelType poly = 0;
  for(size_t i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> poly;
    tetraList[i].setPoly(poly);
  }

  meshFile.close();
}

void MeshReaderPoly::setSections(std::array<std::string, 4> sections)
{
  sections_ = sections;
}

std::array<std::string, 4> MeshReaderPoly::getSections() const
{
  return sections_;
}

void MeshReaderPoly::goToSection(std::ifstream& meshFile, unsigned secNo) const
{
  std::string curLine = "";
  while(curLine != sections_[secNo])
    meshFile >> curLine;
}

}
