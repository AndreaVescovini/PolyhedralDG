#include "MeshReaderPoly.hpp"

#include <vector>
#include <fstream>
#include "Point.hpp"
#include "Tetrahedron.hpp"

namespace dgfem {

MeshReaderPoly::MeshReaderPoly(std::string folder, std::string titleSec1,
                       std::string titleSec2, std::string titleSec3)
  : MeshReader(folder), sections_{{titleSec1, titleSec2, titleSec3}} {}

void MeshReaderPoly::read(Mesh& mesh, const std::string& fileName) const
{
  std::ifstream meshFile{folder_ + fileName};
  if( meshFile.is_open() == false )
  {
    std::cerr << "Can't open mesh file. Mesh file does not exist or is corrupted"
              << std::endl;
  }

  // Use the proxy to access the Mesh class
  MeshProxy mp(mesh);
  std::vector<geom::Point>& vertList = mp.getVerticesRef();
  std::vector<geom::Tetrahedron>& tetraList = mp.getTetrahedraRef();

  // Skip the introduction
  std::string curLine;
  while( curLine != sections_[0] )
    meshFile >> curLine;

  // Read vertices
  size_t verticesNo = 0;
  meshFile >> verticesNo;
  vertList.reserve(verticesNo);

  std::array<real, 3> curVertex;
  unsigned tmp = 0;
  for(size_t i = 0; i < verticesNo; i++)
  {
    meshFile >> curVertex[0];
    meshFile >> curVertex[1];
    meshFile >> curVertex[2];
    vertList.emplace_back(geom::Point(curVertex));
    meshFile >> tmp; // label of the vertex
  }

  // Skip middle lines
  while( curLine != sections_[1] )
    meshFile >> curLine;

  // Read tetrahedra
  size_t tetrahedraNo = 0;
  meshFile >> tetrahedraNo;
  tetraList.reserve(tetrahedraNo);

  std::array<labelType, 4> curTetra;
  for(size_t i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> curTetra[0];
    meshFile >> curTetra[1];
    meshFile >> curTetra[2];
    meshFile >> curTetra[3];
    tetraList.emplace_back(geom::Tetrahedron(curTetra));
    meshFile >> tmp; // label of the vertex
  }

  // Skip middle lines
  while( curLine != sections_[2] )
    meshFile >> curLine;

  // Read polyhedra
  size_t polyhedraNo = 0;
  meshFile >> polyhedraNo;
  labelType poly = 0;
  for(size_t i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> poly;
    tetraList[i].setPoly(poly);
  }

  meshFile.close();
}

void MeshReaderPoly::setSections(std::string titleSec1, std::string titleSec2,
                                 std::string titleSec3)
{
  sections_ = {{titleSec1, titleSec2, titleSec3}};
}

std::array<std::string, 3> MeshReaderPoly::getSections() const
{
  return sections_;
}

}
