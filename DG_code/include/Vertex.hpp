#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include "PolyDG.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <functional>
#include <iostream>

namespace PolyDG
{

class Vertex
{
public:
  Vertex(Real x = 0.0, Real y = 0.0, Real z = 0.0);

  // Default copy constructor.
  Vertex(const Vertex&) = default;

  // Default move constructor.
  Vertex(Vertex&&) = default;

  inline const Eigen::Vector3d& getCoords() const;

  template <typename D>
  void setCoords(const Eigen::MatrixBase<D>& coords);

  inline Real getX() const;
  inline Real getY() const;
  inline Real getZ() const;
  inline unsigned getId() const;

  // Function that returns the distance with v2
  inline Real distance(const Vertex& v2) const;

  // Function that resets the counter (to be used at the beginning of the
  // reading of a new mesh, to assure that ids start from 0).
  inline static void resetCounter(unsigned counter = 0);

  virtual ~Vertex() = default;

  friend std::ostream& operator<<(std::ostream& out, const Vertex& v);

  // Binary comparison operators used to sort face vertices.
  inline friend bool compId(const Vertex& lhs, const Vertex& rhs);
  inline friend bool operator==(const Vertex& lhs, const Vertex& rhs);

private:
  const unsigned id_;
  Eigen::Vector3d coords_;

  // Counter used to assign a different id to each vertex when it is created.
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Vertex& v);

} // namespace PolyDG

namespace std
{

// I specify the equal_to and hash structs in order to use in the right way
// unordered sets of vertices, comparing only the index to see weather two
// elements are equivalent.
template<>
struct equal_to<PolyDG::Vertex>
{
  bool operator()(const PolyDG::Vertex& lhs, const PolyDG::Vertex& rhs) const
  {
    return lhs == rhs;
  }
};

template<>
struct hash<PolyDG::Vertex>
{
  size_t operator()(const PolyDG::Vertex& v) const
  {
    return hash<unsigned>()(v.getId());
  }
};

} // namespace std

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

namespace PolyDG
{

inline const Eigen::Vector3d& Vertex::getCoords() const
{
  return coords_;
}

inline Real Vertex::getX() const
{
  return coords_(0);
}

inline Real Vertex::getY() const
{
  return coords_(1);
}

inline Real Vertex::getZ() const
{
  return coords_(2);
}

inline unsigned Vertex::getId() const
{
  return id_;
}

inline Real Vertex::distance(const Vertex& v2) const
{
  return (coords_ - v2.getCoords()).norm();
}

inline void Vertex::resetCounter(unsigned counter)
{
  counter_ = counter;
}

template <typename D>
void Vertex::setCoords(const Eigen::MatrixBase<D>& coords)
{
  coords_ = coords;
}

inline bool compId(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.id_ < rhs.id_;
}

inline bool operator==(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.id_ == rhs.id_;
}

} // namespace PolyDG

#endif // _VERTEX_HPP_
