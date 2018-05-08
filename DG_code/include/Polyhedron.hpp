#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include "PolyDG.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <Eigen/Geometry>

#include <array>
#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

namespace PolyDG
{

class Tetrahedron;

class Polyhedron
{
public:
  // The defolut constructor also sets the id using the counter.
  Polyhedron();

  // Default copy constructor.
  Polyhedron(const Polyhedron&) = default;

  // Default move constructor.
  Polyhedron(Polyhedron&&) = default;

  inline void addTetra(Tetrahedron& tet);
  inline void addVertexExt(Vertex& v);

  inline unsigned getId() const;
  inline const Tetrahedron& getTetra(SizeType i) const;
  inline Tetrahedron& getTetra(SizeType i);
  inline SizeType getTetrahedraNo() const;
  inline SizeType getVerticesExtNo() const;
  inline const Eigen::AlignedBox3d& getBoundingBox() const;
  inline Real getDiameter() const;

  // Function that starting from the unordered set of vertices computes the
  // cartesian bounding box of the polyhedron.
  void computeBB();

  // Function that starting from the unordered set of vertices computes the
  // diameter of the polyhderon i.e. the maximum distance between two points
  // inside the polyhderon.
  void computeDiameter();

  inline static void resetCounter(unsigned counter = 0);

  virtual ~Polyhedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Polyhedron& poly);

private:
  // Id number of the polyhedron.
  const unsigned id_;

  // Vector containing the tetrahedron of which the polyhedron is made.
  std::vector<std::reference_wrapper<Tetrahedron>> tetrahedra_;

  // Unordered set containing the vertices coming from faces of the polyhedron.
  std::unordered_set<std::reference_wrapper<Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> verticesExt_;

  // Cartesian bounding box containing the polyhedron.
  Eigen::AlignedBox3d boundingBox_;

  // Diameter of the polyhedron i.e. maximum distance between two vertices.
  Real diameter_;

  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Polyhedron& poly);


//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline void Polyhedron::addTetra(Tetrahedron& tet)
{
  tetrahedra_.emplace_back(tet);
}

inline void Polyhedron::addVertexExt(Vertex& v)
{
  verticesExt_.emplace(v);
}

inline Real Polyhedron::getDiameter() const
{
  return diameter_;
}

inline const Tetrahedron& Polyhedron::getTetra(SizeType i) const
{
  return tetrahedra_[i];
}

inline Tetrahedron& Polyhedron::getTetra(SizeType i)
{
  return tetrahedra_[i];
}

inline SizeType Polyhedron::getTetrahedraNo() const
{
  return tetrahedra_.size();
}

inline unsigned Polyhedron::getId() const
{
  return id_;
}

inline const Eigen::AlignedBox3d& Polyhedron::getBoundingBox() const
{
  return boundingBox_;
}

inline SizeType Polyhedron::getVerticesExtNo() const
{
  return verticesExt_.size();
}

inline void Polyhedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _POLYHEDRON_HPP_
