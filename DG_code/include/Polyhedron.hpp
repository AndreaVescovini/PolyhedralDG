#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include "PolyDG.hpp"
#include "Vertex.hpp"
#include "Tetrahedron.hpp"

#include <Eigen/Geometry>

#include <vector>
#include <array>
#include <functional>
#include <unordered_set>

namespace PolyDG
{

class Tetrahedron;

class Polyhedron
{
public:
  // The defolut constructor also sets the id using the counter.
  Polyhedron();

  inline void addTetra(Tetrahedron& tet);
  inline void addVertex(Vertex& v);

  inline Real getDiameter() const;
  inline const Tetrahedron& getTetra(unsigned i) const;
  inline Tetrahedron& getTetra(unsigned i);
  inline unsigned getTetrahedraNo() const;
  inline unsigned getId() const;
  inline const Eigen::AlignedBox3d& getBoundingBox() const;

  // Function that starting from the unordered set of vertices computes the
  // cartesian bounding box of the polyhedron.
  void computeBB();

  // Function that starting from the unordered set of vertices computes the
  // diameter of the polyhderon i.e. the maximum distance between two points
  // inside the polyhderon.
  void computeDiameter();

  inline static void resetCounter(unsigned counter = 0);

  virtual ~Polyhedron() = default;

private:
  const unsigned id_;
  std::vector<std::reference_wrapper<Tetrahedron>> tetrahedra_;

  // Here I store vertices coming from faces
  std::unordered_set<std::reference_wrapper<Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> vertices_;
  // std::array<interval, 3> boundingBox_;
  Eigen::AlignedBox3d boundingBox_;
  Real diameter_;

  static unsigned counter_;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline void Polyhedron::addTetra(Tetrahedron& tet)
{
  tetrahedra_.emplace_back(tet);
}

inline void Polyhedron::addVertex(Vertex& v)
{
  vertices_.emplace(v);
}

inline Real Polyhedron::getDiameter() const
{
  return diameter_;
}

inline const Tetrahedron& Polyhedron::getTetra(unsigned i) const
{
  return tetrahedra_[i];
}

inline Tetrahedron& Polyhedron::getTetra(unsigned i)
{
  return tetrahedra_[i];
}

inline unsigned Polyhedron::getTetrahedraNo() const
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

inline void Polyhedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _POLYHEDRON_HPP_
