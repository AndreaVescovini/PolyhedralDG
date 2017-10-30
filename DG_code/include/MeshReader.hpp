#ifndef _MESH_READER_HPP_
#define _MESH_READER_HPP_

#include <string>
#include "Mesh.hpp"

class MeshReader
{
public:
  MeshReader(std::string folder, std::string titleSec1,
             std::string titleSec2, std::string titleSec3);

  void read(geom::Mesh& Th, const std::string& fileName) const;

private:
  std::string folder_;
  std::string titleSec1_;
  std::string titleSec2_;
  std::string titleSec3_;
};

#endif // _MESH_READER_HPP_
