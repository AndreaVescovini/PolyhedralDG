#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "geom.hpp"

namespace geom {

class Vertex
{
public:
  // Vertex();
  Vertex(real x, real y, real z);

  Vertex(const Vertex&) = default;
  Vertex& operator=(const Vertex&) = default;
  Vertex(Vertex&&) = default;
  Vertex& operator=(Vertex&&) = default;

  const Eigen::Vector3d& getCoords() const;
  real getX() const;
  real getY() const;
  real getZ() const;
  unsigned getId() const;

  real distance(const Vertex& v2) const;

  static void resetCounter(unsigned counter = 0);

  virtual ~Vertex() = default;

  friend std::ostream& operator<<(std::ostream& out, const Vertex& v);
  friend bool compX(const Vertex& lhs, const Vertex& rhs);
  friend bool compY(const Vertex& lhs, const Vertex& rhs);
  friend bool compZ(const Vertex& lhs, const Vertex& rhs);

private:
  const unsigned id_;
  // const std::array<real, 3> coords_;
  const Eigen::Vector3d coords_;

  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Vertex& v);
bool compX(const Vertex& lhs, const Vertex& rhs);
bool compY(const Vertex& lhs, const Vertex& rhs);
bool compZ(const Vertex& lhs, const Vertex& rhs);

}

namespace std {

using geom::Vertex;

template<>
struct equal_to<Vertex>
{
  bool operator()(const Vertex& lhs, const Vertex& rhs) const
  {
    return equal_to<unsigned>()(lhs.getId(), rhs.getId());
  }
};

template<>
struct hash<Vertex>
{
  std::size_t operator()(const Vertex& v) const
  {
    return hash<unsigned>()(v.getId());
  }
};

}

#endif // _VERTEX_HPP_
