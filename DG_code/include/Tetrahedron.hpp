#ifndef _TETRAHEDRON_HPP_
#define _TETRAHEDRON_HPP_

#include <array>
#include <iostream>
#include <Eigen/Geometry>
#include "geom.hpp"
#include "Vertex.hpp"
#include "Polyhedron.hpp"

namespace geom {

class Polyhedron;

class Tetrahedron
{
public:
  Tetrahedron(Vertex& v1,  Vertex& v2, Vertex& v3, Vertex& v4,
              Polyhedron* poly = nullptr);
  Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4,
              Polyhedron& poly);

  Tetrahedron(const Tetrahedron&) = default;
  Tetrahedron& operator=(const Tetrahedron&) = default;
  Tetrahedron(Tetrahedron&&) = default;
  Tetrahedron& operator=(Tetrahedron&&) = default;

  inline const Vertex& getVertex(unsigned i) const;
  inline Vertex& getVertex(unsigned i);
  inline const Polyhedron& getPoly() const;// forse questo metodo è da togliere
  inline Polyhedron& getPoly(); // forse questo const è da togliere
  inline void setPoly(Polyhedron* poly);
  inline void setPoly(Polyhedron& poly);
  inline unsigned getId() const;

  inline const Eigen::Transform<real, 3, Eigen::AffineCompact>& getMap() const;
  inline real getAbsDetJacobian() const;

  static void resetCounter(unsigned counter = 0);

  virtual ~Tetrahedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  const unsigned id_;
  std::array<std::reference_wrapper<Vertex>, 4> vertices_;
  Polyhedron* poly_;
  Eigen::Transform<real, 3, Eigen::AffineCompact> map_;
  real absDetJacobian_;
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

inline void Tetrahedron::setPoly(Polyhedron* poly)
{
  poly_ = poly;
}

inline void Tetrahedron::setPoly(Polyhedron& poly)
{
  this->setPoly(&poly);
}

inline unsigned Tetrahedron::getId() const
{
  return id_;
}

inline const Eigen::Transform<real, 3, Eigen::AffineCompact>& Tetrahedron::getMap() const
{
  return map_;
}

inline real Tetrahedron::getAbsDetJacobian() const
{
  return absDetJacobian_;
}

} // namespace geom

#endif // _TETRAHEDRON_HPP_
