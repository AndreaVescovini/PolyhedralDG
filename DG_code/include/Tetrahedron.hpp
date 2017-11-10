#ifndef _TETRAHEDRON_HPP_
#define _TETRAHEDRON_HPP_

#include <array>
#include <iostream>
#include "geom.hpp"
#include "Vertex.hpp"
#include "Polyhedron.hpp"

namespace geom {

class Polyhedron;

class Tetrahedron
{
public:
  // Tetrahedron() = default;
  Tetrahedron(const Vertex& v1, const Vertex& v2, const Vertex& v3,
              const Vertex& v4, const Polyhedron* poly = nullptr);
  Tetrahedron(const Vertex& v1, const Vertex& v2, const Vertex& v3,
              const Vertex& v4, const Polyhedron& poly);

  Tetrahedron(const Tetrahedron&) = default;
  Tetrahedron& operator=(const Tetrahedron&) = default;
  Tetrahedron(Tetrahedron&&) = default;
  Tetrahedron& operator=(Tetrahedron&&) = default;

  const Vertex& getVertex(unsigned i) const;
  const Polyhedron& getPoly() const;
  void setPoly(const Polyhedron* poly);
  void setPoly(const Polyhedron& poly);
  unsigned getId() const;

  static void resetCounter(unsigned counter = 0);

  virtual ~Tetrahedron() = default;

  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  const unsigned id_;
  const std::array<std::reference_wrapper<const Vertex>, 4> vertices_;
  Polyhedron const* poly_;

  static unsigned counter_;
};

}

#endif // _TETRAHEDRON_HPP_
