#include "MeshReader.hpp"

namespace dgfem {

MeshReader::MeshReader(std::string folder)
  : folder_{folder} {}

void MeshReader::setFolder(const std::string& folder)
{
  folder_ = folder;
}

std::string MeshReader::getFolder() const
{
  return folder_;
}

}
