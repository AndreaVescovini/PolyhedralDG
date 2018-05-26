/*!
    @file   Vertex.hpp
    @author Andrea Vescovini
    @brief  Class that defines a vertex of a tridimensional mesh
*/

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

/*!
    @brief Class that defines a vertex of a tridimensional mesh

    This class defines a vertex of a tridimensional mesh. It stores the
    coordinates and a id number univocal inside the mesh.@n
    A method Vertex::distance is provided for computing the euclidean distance
    between two vertices.
*/

class Vertex
{
public:
  /*!
      @brief Constructor that takes the three coordinates

      This constructor takes three Reals as the three coordinates. If left empty
      the coordinates are set to zero.
  */
  Vertex(Real x = 0.0, Real y = 0.0, Real z = 0.0);

  //! Copy constructor
  Vertex(const Vertex&) = default;

  //! Move constructor
  Vertex(Vertex&&) = default;

  /*!
      @brief  Get the coordinates
      @return A const reference to the coordinates.
  */
  inline const Eigen::Vector3d& getCoords() const;

  /*!
      @brief Set the coordinates
      @param coords Vector containing the coordinates.
  */
  template <typename D>
  void setCoords(const Eigen::MatrixBase<D>& coords);

  //! Get the x coordinate
  inline Real getX() const;

  //! Get the y coordinate
  inline Real getY() const;

  //! Get the z coordinate
  inline Real getZ() const;

  //! Get the id number
  inline unsigned getId() const;

  //! Compute the euclidean distance with another Vertex
  inline Real distance(const Vertex& v2) const;

  /*!
      @brief Reset the counter for id numbers

      This function resets the counter that assigns id numbers. It should be used
      before reading of a new mesh, to assure that id numbers start from 0.

      @param counter The number to which the counter is reset. If not specified
                     it is 0.
  */
  inline static void resetCounter(unsigned counter = 0);

  //! Destructor
  virtual ~Vertex() = default;

  /*!
      @brief Overload for the ostream operator

      It prints the id number followed by the coordinates.
  */
  friend std::ostream& operator<<(std::ostream& out, const Vertex& v);

  //! Comparison operator< based on the id number
  inline friend bool compareId(const Vertex& lhs, const Vertex& rhs);

  //! Overload for the operator== based on the id number
  inline friend bool operator==(const Vertex& lhs, const Vertex& rhs);

private:
  //! Id number
  const unsigned id_;

  //! Coordinates
  Eigen::Vector3d coords_;

  //! Counter used to assign a different id to each vertex when it is created
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Vertex& v);

} // namespace PolyDG

namespace std
{

// I specify the equal_to and hash structs in order to use in the right way
// unordered sets of vertices, comparing only the index to see weather two
// elements are equivalent.
//! Specialization of the functor @c std::equal_to for a PolyDG::Vertex
template<>
struct equal_to<PolyDG::Vertex>
{
  //! Call operator
  bool operator()(const PolyDG::Vertex& lhs, const PolyDG::Vertex& rhs) const
  {
    return lhs == rhs;
  }
};

//! Specialization of the functor @c std::hash for a PolyDG::Vertex based of the id number of the Vertex
template<>
struct hash<PolyDG::Vertex>
{
  //! Call operator
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

//! Comparison operator< based on the id number
inline bool compareId(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.id_ < rhs.id_;
}

//! Overload for the operator== based on the id number
inline bool operator==(const Vertex& lhs, const Vertex& rhs)
{
  return lhs.id_ == rhs.id_;
}

} // namespace PolyDG

#endif // _VERTEX_HPP_
