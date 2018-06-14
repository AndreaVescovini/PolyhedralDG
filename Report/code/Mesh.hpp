class Mesh
{
public:
  Mesh(const std::string& fileName, MeshReader& reader)
  {
    reader.read(*this, fileName);
    // ...
    computeFaces();
    // ...
  }

  // Getter methods
  inline const Vertex& getVertex(SizeType i) const;
  inline SizeType getVerticesNo() const;
  // The same for the other entities...

  // ...

  // Proxy class used to modify the mesh
  friend class MeshProxy;

private:
  std::vector<Vertex> vertices_;
  std::vector<Tetrahedron> tetrahedra_;
  std::vector<FaceExt> facesExt_;
  std::vector<FaceInt> facesInt_;
  std::vector<Polyhedron> polyhedra_;

  // Compute the internal faces of the mesh and complete the information about the external ones
  void computeFaces();

  // ...
};
