/*!
    @file   Tetrahedron.hpp
    @author Andrea Vescovini
    @brief  Class that defines a tetrahedron of a tridimensional mesh
*/

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

/*!
    @brief Class that defines a tetrahedron of a tridimensional mesh

    This class defines a tetrahedron of a tridimensional mesh. It stores the
    four vertices, a id number univocal inside the mesh, the affine map from the
    reference tetrahedron (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1) to this and
    the absolute value of the determinant of the jacobian of the linear part of
    the map. It stores also a pointer to a Polyhedron to which the tetrahedron belongs.
*/

class Tetrahedron
{
public:
  /*!
      @brief Constructor that takes four vertices

      The four vertices are setted and the map. The absolute value of the
      determinant of its jacobian are initializated.
  */
  Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4);

  /*!
      @brief Constructor that takes four vertices and a reference to a Polyhedron

      The four vertices and the Polyhedron are setted. The map and the absolute
      value of the determinant of its jacobian are initializated.
  */
  Tetrahedron(Vertex& v1, Vertex& v2, Vertex& v3, Vertex& v4, Polyhedron& poly);

  //! Copy constructor
  Tetrahedron(const Tetrahedron&) = default;

  //! Move constructor
  Tetrahedron(Tetrahedron&&) = default;

  /*!
      @brief Get a Vertex

      This functions returns the i-th Vertex.

      @param i The index of the Vertex required, it can be 0, 1, 2 or 3.
  */
  inline const Vertex& getVertex(unsigned i) const;

  /*!
      @brief Get a Vertex

      This functions returns the i-th Vertex.

      @param i The index of the Vertex required, it can be 0, 1, 2 or 3.
  */
  inline Vertex& getVertex(unsigned i);

  /*!
      @brief  Check if a Polyhedron is set
      @return @b true if a Polyhedron is set, @b false if it is not.
  */
  inline bool isPolySet() const;

  /*!
      @brief   Get the Polyhedron to which the Tetrahedron belongs
      @warning This function must be called only if isPolySet() returns @b true.
  */
  inline const Polyhedron& getPoly() const;

  /*!
      @brief   Get the Polyhedron to which the Tetrahedron belongs
      @warning This function must be called only if isPolySet() returns @b true.
  */
  inline Polyhedron& getPoly();

  /*!
      @brief Set a Polyhedron
      @param poly The Polyhedron to which this Tetrahedron belongs.
  */
  inline void setPoly(Polyhedron& poly);

  //! Get the id number
  inline unsigned getId() const;

  /*!
      @brief  Get the affine map from the reference tetrahedron
      @return a Eigen::Transform object containing the affine map from the
              reference tetrahedron (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
              to this Tetrahedron.
  */
  inline const Eigen::Transform<Real, 3, Eigen::AffineCompact>& getMap() const;

  /*!
      @brief  Get the absolute value of the determinant of the jacobian
      @return The absolute value of the determinant of the jacobian of the
              affine map from the reference tetrahedron (0, 0, 0), (1, 0, 0),
              (0, 1, 0), (0, 0, 1) to this Tetrahedron.
  */
  inline Real getAbsDetJacobian() const;

  /*!
      @brief Reset the counter for id numbers

      This function resets the counter that assigns id numbers. It should be used
      before reading of a new mesh, to assure that id numbers start from 0.

      @param counter The number to which the counter is reset. If not specified
                     it is 0.
  */
  inline static void resetCounter(unsigned counter = 0);

  //! Destructor
  virtual ~Tetrahedron() = default;

  /*!
      @brief Overload for the operator<<

      It prints the id number followed by the id number of the four vertices
      and the id number of the Polyhedron, if it is set.
  */
  friend std::ostream& operator<<(std::ostream& out, const Tetrahedron& tetra);

private:
  //! id number
  const unsigned id_;

  //! The four vertices of the tetrahedron
  std::array<std::reference_wrapper<Vertex>, 4> vertices_;

  //! The Polyhedron to which the tetrahedron belongs
  Polyhedron* poly_;

  //! The affine map from the reference tetrahedron
  Eigen::Transform<Real, 3, Eigen::AffineCompact> map_;

  //! The absolute value of the determinant of the jacobian of the map
  Real absDetJacobian_;

  //! Counter used to assign a different id to each vertex when it is created.
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

inline bool Tetrahedron::isPolySet() const
{
  return poly_ == nullptr ? false : true;
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
