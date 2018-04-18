#ifndef _MESH_READER_HPP_
#define _MESH_READER_HPP_

#include "Mesh.hpp"

#include <string>

namespace PolyDG
{

class Mesh;

class MeshReader
{
public:
  MeshReader() = default;

  virtual void read(Mesh& mesh, const std::string& fileName) const = 0;

  virtual ~MeshReader() = default;
};

} // namespace PolyDG

#endif // _MESH_READER_HPP_
