/*!
    @file   Mesh.hpp
    @author Andrea Vescovini
    @brief  Class that defines a polyhedral mesh
*/

#ifndef _MESH_HPP_
#define _MESH_HPP_

#include "FaceExt.hpp"
#include "FaceInt.hpp"
#include "MeshProxy.hpp"
#include "MeshReader.hpp"
#include "PolyDG.hpp"
#include "Polyhedron.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PolyDG
{

class MeshReader;

/*!
    @brief Class that defines a polyhedral mesh

    This class defines a polyhedral mesh. It stores vertices, tetrahedra,
    polyhedra, external faces and internal faces. A Mesh can be constructed
    passing to the constructor the name of the file to be read and a MeshReader
    that, through a MeshProxy, fills the containers of the Mesh.
*/
class Mesh
{
public:
  /*!
      @brief Constructor

      This constructor intializes the Mesh reading the file fileName through the
      MeshReader reader, then computes the internal faces of the Mesh if they
      have not been provided and at last the maximum and minimum diameter.

      @param fileName std::string containing the name of the file with the mesh
                      to be read.
      @param reader   MeshReader that provides a function MeshReader::read() which
                      interprets the file with the mesh and stores it in the Mesh.

  */
  Mesh(const std::string& fileName, MeshReader& reader);

  //! Copy constructor
  Mesh(const Mesh&) = default;

  //! Copy-assigment operator
  Mesh& operator=(const Mesh&) = default;

  //! Move constructor
  Mesh(Mesh&&) = default;

  //! Move-assignment operator
  Mesh& operator=(Mesh&&) = default;

  /*!
      @brief Get a Vertex

      This functions returns the i-th Vertex.

      @param i The index of the Vertex required, it can be 0,..,geVerticesNo() - 1.
  */
  inline const Vertex& getVertex(SizeType i) const;

  /*!
      @brief Get a Tetrahedron

      This functions returns the i-th Tetrahedron.

      @param i The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1.
  */
  inline const Tetrahedron& getTetrahedron(SizeType i) const;

  /*!
      @brief Get a FaceExt

      This functions returns the i-th external face.

      @param i The index of the FaceExt required, it can be 0,..,getFacesExtNo() - 1.
  */
  inline const FaceExt& getFaceExt(SizeType i) const;

  /*!
      @brief Get a FaceInt

      This functions returns the i-th internal face.

      @param i The index of the FaceInt required, it can be 0,..,getFacesIntNo() - 1.
  */
  inline const FaceInt& getFaceInt(SizeType i) const;

  /*!
      @brief Get a Polyhedron

      This functions returns the i-th Polyhedron.

      @param i The index of the Polyhedron required, it can be 0,..,getPolyhedraNo() - 1.
  */
  inline const Polyhedron& getPolyhedron(SizeType i) const;

  //! Get the number of vetices
  inline SizeType getVerticesNo() const;

  //! Get the number of tetrahedra
  inline SizeType getTetrahedraNo() const;

  //! Get the number of external faces
  inline SizeType getFacesExtNo() const;

  //! Get the number of internal faces
  inline SizeType getFacesIntNo() const;

  //! Get the number of polyhedra
  inline SizeType getPolyhedraNo() const;

  //! Get the maximum diameter of the polyhedra
  inline Real getMaxDiameter() const;

  //! Get the minumum diameter of the polyhedra
  inline Real getMinDiameter() const;

  //! Prints all the entries of each entity of the Mesh
  void printAll(std::ostream& out = std::cout) const;

  //! Prints the first five entries of each entity of the Mesh
  void printHead(std::ostream& out = std::cout) const;

  //! Prints general information about the Mesh
  void printInfo(std::ostream& out = std::cout) const;

  /*!
      @brief Export the Mesh

      This function exports the mesh into a VTK unstructured grid file with
      XML format. It can be read with a visualization software (e.g. Paraview).

      @param fileName  Name of the file to be saved (the extension should be .vtu).
      @param precision Precision to be used for floating points numbers.
  */
  void exportMeshVTK(const std::string& fileName, unsigned precision = 8) const;

  //! Destructor
  virtual ~Mesh() = default;

  //! Proxy class used to modify the mesh
  friend class MeshProxy;

private:
  //! Vector of Vertex
  std::vector<Vertex> vertices_;

  //! Vector of Tetrahedron
  std::vector<Tetrahedron> tetrahedra_;

  //! Vector of external faces
  std::vector<FaceExt> facesExt_;

  //! Vector of internal faces
  std::vector<FaceInt> facesInt_;

  //! Vector of polyhedra
  std::vector<Polyhedron> polyhedra_;

  //! Maximum diameter
  Real hmax_;

  //! Minumim diameter
  Real hmin_;

  //! Computes the internal faces of the mesh and complete the information about the external ones
  void computeFaces();

  //! Inialize the maximum and minumum diameter
  void computeDiameters();

  //! Prints the first lineNo entries of each entity of the mesh
  void print(SizeType lineNo, std::ostream& out = std::cout) const;
};

/*!
    @brief Exception to be thrown in case of error in the file with the mesh

    This class inherits from @c std::runtime_error and should be used in case of
    discrepancies beteen the file with the mesh and the behaviour of the
    MeshReader.
*/
class MeshFormatError : public std::runtime_error
{
public:
  //! Constructor that takes a std::string
  explicit MeshFormatError(const std::string& what_arg);

  //! Constructor that takes a C-style string
  explicit MeshFormatError(const char* what_arg);

  //! Destructor
  virtual ~MeshFormatError() = default;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const Vertex& Mesh::getVertex(SizeType i) const
{
  return vertices_[i];
}

inline const Tetrahedron& Mesh::getTetrahedron(SizeType i) const
{
  return tetrahedra_[i];
}

inline const FaceExt& Mesh::getFaceExt(SizeType i) const
{
  return facesExt_[i];
}

inline const FaceInt& Mesh::getFaceInt(SizeType i) const
{
  return facesInt_[i];
}

inline const Polyhedron& Mesh::getPolyhedron(SizeType i) const
{
  return polyhedra_[i];
}

inline SizeType Mesh::getVerticesNo() const
{
  return vertices_.size();
}

inline SizeType Mesh::getTetrahedraNo() const
{
  return tetrahedra_.size();
}

inline SizeType Mesh::getFacesExtNo() const
{
  return facesExt_.size();
}

inline SizeType Mesh::getFacesIntNo() const
{
  return facesInt_.size();
}

inline SizeType Mesh::getPolyhedraNo() const
{
  return polyhedra_.size();
}

inline Real Mesh::getMaxDiameter() const
{
  return hmax_;
}

inline Real Mesh::getMinDiameter() const
{
  return hmin_;
}

} // namespace PolyDG

#endif // _MESH_HPP_
