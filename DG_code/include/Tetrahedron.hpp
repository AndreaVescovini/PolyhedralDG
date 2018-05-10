#ifndef _TETRAHEDRON_HPP_
#define _TETRAHEDRON_HPP_

#include "PolyDG.hpp"
#include "Polyhedron.hpp"
#include "Vertex.hpp"

#include <Eigen/Geometry>

#include <array>
#include <iostream>

namespace PolyDG
{

class Polyhedron;

class Tetrahedron
{
public:
  Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4);
  Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4, Polyhedron& poly);

  // Default copy constructor.
  Tetrahedron(const Tetrahedron&) = default;

  // Default move constructor.
  Tetrahedron(Tetrahedron&&) = default;

  inline const Vertex& getVertex(unsigned i) const;
  inline Vertex& getVertex(unsigned i);
  inline const Polyhedron& getPoly() const;
  inline Polyhedron& getPoly();
  inline void setPoly(Polyhedron& poly);
  inline unsigned getId() const;

  inline const Eigen::Transform<Real, 3, Eigen::AffineCompact>& getMap() const;
  inline Real getAbsDetJacobian() const;

  inline static void resetCounter(unsigned counter = 0);

  virtual ~Tetrahedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  const unsigned id_;
  std::array<std::reference_wrapper<Vertex>, 4> vertices_;
  Polyhedron* poly_;
  Eigen::Transform<Real, 3, Eigen::AffineCompact> map_;
  Real absDetJacobian_;
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline const Vertex& Tetrahedron::getVertex(unsigned i) const
{
  return vertices_[i];
}

inline Vertex& Tetrahedron::getVertex(unsigned i)
{
  return vertices_[i];
}

inline const Polyhedron& Tetrahedron::getPoly() const
{
  return *poly_;
}

inline Polyhedron& Tetrahedron::getPoly()
{
  return *poly_;
}

inline void Tetrahedron::setPoly(Polyhedron& poly)
{
  poly_ = &poly;
}

inline unsigned Tetrahedron::getId() const
{
  return id_;
}

inline const Eigen::Transform<Real, 3, Eigen::AffineCompact>& Tetrahedron::getMap() const
{
  return map_;
}

inline Real Tetrahedron::getAbsDetJacobian() const
{
  return absDetJacobian_;
}

inline void Tetrahedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _TETRAHEDRON_HPP_
