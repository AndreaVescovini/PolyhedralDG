#include "MeshReader.hpp"

#include <vector>
#include <array>
#include <fstream>
#include "geom.hpp"
#include "Point.hpp"
#include "Tetrahedron.hpp"

using geom::Point;
using geom::Tetrahedron;
using geom::Mesh;

MeshReader::MeshReader(std::string folder, std::string titleSec1,
                       std::string titleSec2, std::string titleSec3)
  : folder_{folder}, titleSec1_{titleSec1}, titleSec2_{titleSec2},
    titleSec3_{titleSec3} {}

void MeshReader::read(Mesh& Th, const std::string& fileName) const
{
  std::ifstream meshFile{folder_ + fileName};
  if( meshFile.is_open() == false )
  {
    std::cerr << "Can't open mesh file. Mesh file does not exist or is corrupted"
              << std::endl;
  }

  // Skip the introduction
  std::string curLine;
  while( curLine != titleSec1_ )
    meshFile >> curLine;

  // Read vertices
  unsigned verticesNo = 0;
  meshFile >> verticesNo;
  std::vector<Point> vertices(verticesNo);

  std::array<real, 3> curVertex;
  unsigned tmp = 0;
  for(unsigned i = 0; i < verticesNo; i++)
  {
    meshFile >> curVertex[0];
    meshFile >> curVertex[1];
    meshFile >> curVertex[2];
    vertices[i] = Point(curVertex);
    meshFile >> tmp; // label of the vertex
  }

  Th.setVerticesNo(verticesNo);
  Th.setVertices(vertices);

  // Skip middle lines
  while( curLine != titleSec2_ )
    meshFile >> curLine;

  // Read tetrahedra
  unsigned tetrahedraNo = 0;
  meshFile >> tetrahedraNo;
  std::vector<Tetrahedron> tetrahedra(tetrahedraNo);

  std::array<labelType, 4> curTetra;
  for(unsigned i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> curTetra[0];
    meshFile >> curTetra[1];
    meshFile >> curTetra[2];
    meshFile >> curTetra[3];
    tetrahedra[i] = Tetrahedron(curTetra);
    meshFile >> tmp; // label of the vertex
  }

  Th.setTetrahedraNo(tetrahedraNo);

  // Skip middle lines
  while( curLine != titleSec3_ )
    meshFile >> curLine;

  // Read polyhedra
  unsigned polyhedraNo = 0;
  meshFile >> polyhedraNo;
  labelType poly = 0;
  for(unsigned i = 0; i < tetrahedraNo; i++)
  {
    meshFile >> poly;
    tetrahedra[i].setPoly(poly);
  }

  Th.setTetrahedra(tetrahedra);
  Th.setPolyhedraNo(polyhedraNo);

  meshFile.close();
}
