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

  const Vertex& getVertex(unsigned i) const;
  Vertex& getVertex(unsigned i);
  const Polyhedron& getPoly() const;// forse questo metodo è da togliere
  Polyhedron& getPoly(); // forse questo const è da togliere
  void setPoly(Polyhedron* poly);
  void setPoly(Polyhedron& poly);
  unsigned getId() const;

  static void resetCounter(unsigned counter = 0);

  virtual ~Tetrahedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  const unsigned id_;
  std::array<std::reference_wrapper<Vertex>, 4> vertices_;
  Polyhedron* poly_;
  Eigen::Transform<real, 3, Eigen::AffineCompact> map_;
  real detJacobian_;
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

}

#endif // _TETRAHEDRON_HPP_
