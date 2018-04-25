#include "Mesh.hpp"
#include "Face.hpp"

#include <algorithm>
#include <unordered_set>
#include <memory>


namespace PolyDG
{

Mesh::Mesh(const std::string& fileName, MeshReader& reader)
{
  reader.read(*this, fileName);
  std::cout << "The mesh file has been read." << std::endl;
  computeFaces();
  std::cout << "Faces have been computed." << std::endl;
  computePolyInfo();
  std::cout << "Informations about polyhedra have been computed." << std::endl;
}

Real Mesh::getMaxDiameter() const
{
  Real hmax = polyhedra_[0].getDiameter();

  for(unsigned i = 1; i < polyhedra_.size(); i++)
    if(polyhedra_[i].getDiameter() > hmax)
      hmax = polyhedra_[i].getDiameter();

  return hmax;
}

void Mesh::printAll(std::ostream& out) const
{
  unsigned lineNo = std::max(vertices_.size(), tetrahedra_.size());
  this->print(lineNo, out);
}

void Mesh::printHead(std::ostream& out) const
{
  this->print(5, out);
}

void Mesh::computeFaces()
{
  std::unordered_set<std::unique_ptr<Face>> temp;
  temp.reserve(tetrahedra_.size() * 4);
  facesInt_.reserve(tetrahedra_.size() * 2);

  for(Tetrahedron& t : tetrahedra_) // perchè non const?
    for(unsigned faceNo = 0; faceNo < 4; faceNo++)
    { // I store faces in temp, the constructor will sort vertices so that
      // the comparison will be easy, if a face is not present it will be inserted,
      // otherwise we will get res.second == false

      // std::cout << "Tetraedro " << t.getId() << " " << faceNo << std::endl;

      // The i-th face is that one without the (3-i)-th vertex.
      auto res = temp.emplace(new Face(t.getVertex(static_cast<unsigned>(faceNo < 1)),
                                       t.getVertex(static_cast<unsigned>(faceNo < 2) + 1),
                                       t.getVertex(static_cast<unsigned>(faceNo < 3) + 2),
                                       t, 3 - faceNo));
      if(res.second == false)
      {
        // If I find a face of a tetrahedron of a different polyhedron I insert it in
        // facesInt_, then in any case I erase the face from temp
        unsigned elem1 = t.getPoly().getId();
        unsigned elem2 = (*res.first)->getTet1().getPoly().getId();

        if(elem1 != elem2)
        {
          if(elem1 > elem2)
            facesInt_.emplace_back((*res.first)->getVertex(0), (*res.first)->getVertex(1), (*res.first)->getVertex(2),
                                   (*res.first)->getTet1(), (*res.first)->getFaceNoTet1(), t, 3 - faceNo);
          else
            facesInt_.emplace_back((*res.first)->getVertex(0), (*res.first)->getVertex(1), (*res.first)->getVertex(2),
                                    t, 3 - faceNo, (*res.first)->getTet1(), (*res.first)->getFaceNoTet1());
        }
        temp.erase(res.first);
      }
    }

  // Now in temp there are only external faces that have been inserted only once,
  // I loop over external faces and I find in the faces stored in temp the
  // information about the tetrahedron to which it belongs and the local number
  // of the face.
  for(FaceExt& f : facesExt_)
  {
    std::unique_ptr<Face> fext(new Face(f));
    auto got = temp.find(fext);
    // Metterci dei move?
    // f.setNormal((*got)->getNormal());
    f.setTet1((*got)->getTet1());
    f.setFaceNoTet1((*got)->getFaceNoTet1());
    // std::cout << "Now I will check the sign of the normal" << std::endl;
    f.checkNormalSign();
  }

}

void Mesh::computePolyInfo()
{
  // I add in each polyhedron the vertices that are on the faces.
  // I will use them in order to compute the bounding box and the diameter.
  for(FaceInt& f : facesInt_)
    for(unsigned i = 0; i < 3; i++)
    {
      f.getTet1().getPoly().addVertexExt(f.getVertex(i));
      f.getTet2().getPoly().addVertexExt(f.getVertex(i));
    }

  for(FaceExt& f : facesExt_)
    for(unsigned i = 0; i < 3; i++)
      f.getTet1().getPoly().addVertexExt(f.getVertex(i));

  for(Polyhedron& p : polyhedra_)
  {
    p.computeBB();
    p.computeDiameter();
  }
}

void Mesh::print(unsigned lineNo, std::ostream& out) const
{
  out << "-------- MESH --------" << std::endl;

  out << "VERTICES: " << vertices_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, vertices_.size()); i++)
    out << vertices_[i] << std::endl;

  out << "\nTETRAHEDRA: " << tetrahedra_.size() << ", POLYHEDRA: " << polyhedra_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, tetrahedra_.size()); i++)
    out << tetrahedra_[i] << std::endl;

  out << "\nEXTERNAL FACES: " << facesExt_.size() << std::endl;

  for(unsigned i = 0; i < std::min<size_t>(lineNo, facesExt_.size()); i++)
    out << facesExt_[i] << std::endl;

  out << "\nINTERNAL FACES: " << facesInt_.size() << std::endl;

    for(unsigned i = 0; i < std::min<size_t>(lineNo, facesInt_.size()); i++)
      out << facesInt_[i] << std::endl;

  out << "----------------------" << std::endl;
}

} // namespace PolyDG
