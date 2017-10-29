#include <fstream>
#include <algorithm>

#include "Mesh.hpp"


namespace geom {

void Mesh::load(const std::string& fileName)
{
  std::string meshDir = "../meshes/";
  std::ifstream f{meshDir + fileName};
  if( f.is_open() == false )
  {
    std::cerr << "Can't open mesh file. Mesh file does not exist or is corrupted"
              << std::endl;
  }

  // Skip the introduction
  std::string curLine;
  while( curLine != "Vertices" )
    f >> curLine;

  // Read vertices
  f >> verticesNo_;
  vertices_.reserve(verticesNo_);

  std::array<real, 3> curVertex;
  unsigned tmp = 0;
  for(unsigned i = 0; i < verticesNo_; i++)
  {
    f >> curVertex[0];
    f >> curVertex[1];
    f >> curVertex[2];
    vertices_.emplace_back(Point(curVertex));
    f >> tmp; // label of the vertex
  }

  // Skip middle lines
  while( curLine != "Tetrahedra" )
    f >> curLine;

  // Read tetrahedra
  f >> tetrahedraNo_;
  tetrahedra_.reserve(tetrahedraNo_);

  std::array<labelType, 4> curTetra;
  for(unsigned i = 0; i < tetrahedraNo_; i++)
  {
    f >> curTetra[0];
    f >> curTetra[1];
    f >> curTetra[2];
    f >> curTetra[3];
    tetrahedra_.emplace_back(Tetrahedron(curTetra));
    f >> tmp; // label of the vertex
  }

  // Skip middle lines
  while( curLine != "Polyhedra" )
    f >> curLine;

  // Read polyhedra
  f >> polyhedraNo_;
  labelType polyNo = 0;
  for(unsigned i = 0; i < tetrahedraNo_; i++)
  {
    f >> polyNo;
    tetrahedra_[i].setPolyNo(polyNo);
  }

  f.close();
}

void Mesh::print(unsigned lineNo, std::ostream& out) const
{
  out << "-------- MESH --------" << std::endl;

  out << "VERTICES: " << verticesNo_ << std::endl;
  for(unsigned i = 0; i < std::min(lineNo, verticesNo_); i++)
    out << vertices_[i] << std::endl;

  out << "\nTETRAHEDRA: " << tetrahedraNo_ << ", POLYHEDRA: " << polyhedraNo_
      << std::endl;
  for(unsigned i = 0; i < std::min(lineNo, tetrahedraNo_); i++)
    out << tetrahedra_[i] << std::endl;

  out << "----------------------" << std::endl;
}

void Mesh::printAll(std::ostream& out) const
{
  unsigned lineNo = std::max(verticesNo_, tetrahedraNo_);
  this->print(lineNo, out);
}

void Mesh::printHead(std::ostream& out) const
{
  this->print(5, out);
}

}
