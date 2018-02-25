#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include <array>
#include <iostream>
#include <functional>
#include <Eigen/Dense>
#include "geom.hpp"

namespace geom {

class Vertex
{
public:
  Vertex(real x = 0.0, real y = 0.0, real z = 0.0);

  Vertex(const Vertex&) = default;
  Vertex& operator=(const Vertex&) = default;
  Vertex(Vertex&&) = default;
  Vertex& operator=(Vertex&&) = default;

  const Eigen::Vector3d& getCoords() const;
  void setCoords(const Eigen::Vector3d& coords);
  real getX() const;
  real getY() const;
  real getZ() const;
  unsigned getId() const;

  real distance(const Vertex& v2) const;

  // Function that resets the counter (to be used at the beginning of the
  // reading of a new mesh, to assure that ids start from 0).
  static void resetCounter(unsigned counter = 0);

  virtual ~Vertex() = default;

  friend std::ostream& operator<<(std::ostream& out, const Vertex& v);

  // Binary comparison operators used to sort face vertices.
  friend bool compX(const Vertex& lhs, const Vertex& rhs);
  friend bool compY(const Vertex& lhs, const Vertex& rhs);
  friend bool compZ(const Vertex& lhs, const Vertex& rhs);
  friend bool compId(const Vertex& lhs, const Vertex& rhs);

private:
  const unsigned id_;
  Eigen::Vector3d coords_;

  // Counter used to assign a different id to each vertex when it is created.
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Vertex& v);
bool compX(const Vertex& lhs, const Vertex& rhs);
bool compY(const Vertex& lhs, const Vertex& rhs);
bool compZ(const Vertex& lhs, const Vertex& rhs);
bool compId(const Vertex& lhs, const Vertex& rhs);

}

namespace std {

using geom::Vertex;

// I specify the equal_to and hash structs in order to use in the right way
// unordered sets of vertices, comparing only the index to see weather two
// elements are equivalent.
template<>
struct equal_to<geom::Vertex>
{
  bool operator()(const geom::Vertex& lhs, const geom::Vertex& rhs) const
  {
    return equal_to<unsigned>()(lhs.getId(), rhs.getId());
  }
};

template<>
struct hash<geom::Vertex>
{
  std::size_t operator()(const geom::Vertex& v) const
  {
    return hash<unsigned>()(v.getId());
  }
};

}

#endif // _VERTEX_HPP_
