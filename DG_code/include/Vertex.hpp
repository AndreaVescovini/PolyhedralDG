#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include <array>
#include <iostream>
#include "geom.hpp"

namespace geom {

class Vertex
{
public:
  // Vertex();
  explicit Vertex(const std::array<real, 3>& coords);

  Vertex(const Vertex&);
  Vertex& operator=(const Vertex&) = default;
  Vertex(Vertex&&) = default;
  Vertex& operator=(Vertex&&) = default;

  std::array<real, 3> getCoords() const;
  real getX() const;
  real getY() const;
  real getZ() const;
  unsigned getId() const;

  real distance(const Vertex& v2) const;

  static void resetCounter(unsigned counter = 0);

  virtual ~Vertex() = default;

  friend std::ostream& operator<<(std::ostream& out, const Vertex& v);

private:
  const unsigned id_;
  const std::array<real, 3> coords_;

  static unsigned counter_;
};

}

#endif // _VERTEX_HPP_
