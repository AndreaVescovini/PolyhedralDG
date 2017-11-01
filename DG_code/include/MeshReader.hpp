#ifndef _MESH_READER_HPP_
#define _MESH_READER_HPP_

#include <string>
#include "Mesh.hpp"

namespace dgfem {

class Mesh;

class MeshReader
{
public:
  MeshReader() = default;
  MeshReader(std::string folder);

  virtual void read(Mesh& mesh, const std::string& fileName) const = 0;

  void setFolder(const std::string& folder);
  std::string getFolder() const;

  virtual ~MeshReader() = default;

protected:
  std::string folder_;
};

}

#endif // _MESH_READER_HPP_
