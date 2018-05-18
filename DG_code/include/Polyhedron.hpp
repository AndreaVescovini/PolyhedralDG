/*!
    @file   Polyhedron.hpp
    @author Andrea Vescovini
    @brief  Class that defines a polyhedron of a polyhedral mesh
*/

#ifndef _POLYHEDRON_HPP_
#define _POLYHEDRON_HPP_

#include "PolyDG.hpp"
#include "Tetrahedron.hpp"
#include "Vertex.hpp"

#include <Eigen/Geometry>

#include <array>
#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

namespace PolyDG
{

class Tetrahedron;

/*!
    @brief Class that defines a polyhedron of a polyhedral mesh

    This class defines a polyhedron that is an element of a polyhedral mesh.
    A polyhedron is defined as the union of disjoint neighbouring tetrahedra.
    Each Polyhedron has a id number univocal inside the mesh.
    The default constructor creates an empty Polyhedron, then through the
    function addTetra() you can add the tetrahedra.
*/

class Polyhedron
{
public:
  //! Alias for a forward const iterator over the vertices of the Polyhedron
  using ConstIterVertex = std::unordered_set<std::reference_wrapper<Vertex>,
                                             std::hash<Vertex>,
                                             std::equal_to<Vertex>>::const_iterator;
  /*!
      @brief Constructor

      This constructor creates an empty Polyhedron and sets the id number.
  */
  Polyhedron();

  //! Copy constructor
  Polyhedron(const Polyhedron&) = default;

  //! Move constructor
  Polyhedron(Polyhedron&&) = default;

  //! Extend the Polyhedron adding a Tetrahedron
  void addTetra(Tetrahedron& tet);

  //! Get the id number
  inline unsigned getId() const;

  /*!
      @brief Get a Tetrahedron

      This functions returns the i-th Tetrahedron.

      @param i The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1
  */
  inline const Tetrahedron& getTetra(SizeType i) const;

  /*!
      @brief Get a Tetrahedron

      This functions returns the i-th Tetrahedron.

      @param i The index of the Tetrahedron required, it can be 0,..,getTetrahedraNo() - 1
  */
  inline Tetrahedron& getTetra(SizeType i);

  //! Get the number of tetrahedra of which the polyhedron is made
  inline SizeType getTetrahedraNo() const;

  /*!
      @brief     Get a ConstIterVertex poining to the first Vertex
      @attention When you dereference the iterator to access the Vertex, you may
                 have to call the method get(), because vertices as stored
                 through std::reference_wrapper.
  */
  inline ConstIterVertex verticesCbegin() const;

  //! Get a ConstIterVertex poining to the @a past-the-end Vertex
  inline ConstIterVertex verticesCend() const;

  //! Get the number of vertices of the tetrahedron
  inline SizeType getVerticesNo() const;

  /*!
      @brief Get the cartesian bounding box
      @return A const reference to an object of type Eigen::AlignedBox3d that
              rapresents the cartesian bounding box containing the Polyhedron.
  */
  inline const Eigen::AlignedBox3d& getBoundingBox() const;

  //! Get the diameter of the Polyhedron
  inline Real getDiameter() const;

  /*!
      @brief Reset the counter for id numbers

      This function resets the counter that assigns id numbers. It should be used
      before reading of a new mesh, to assure that id numbers start from 0.

      @param counter The number to which the counter is reset. If not specified
                     it is 0.
  */
  inline static void resetCounter(unsigned counter = 0);

  //! Destructor
  virtual ~Polyhedron() = default;

  /*!
      @brief Overload for the operator<<

      It prints the id number of the Polyhedron and the the id number of each
      Tetrahedron inside.
  */
  friend std::ostream& operator<<(std::ostream& out, const Polyhedron& poly);

private:
  //! Id number
  const unsigned id_;

  //! Vector containing the tetrahedra of which the polyhedron is made
  std::vector<std::reference_wrapper<Tetrahedron>> tetrahedra_;

  //! Unordered set containing the vertices coming from the tetrahedra
  std::unordered_set<std::reference_wrapper<Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> vertices_;

  //! Cartesian bounding box containing the polyhedron
  Eigen::AlignedBox3d boundingBox_;

  //! Diameter of the polyhedron i.e. maximum distance between two vertices
  Real diameter_;

  //! Counter used to assign a different id to each Polyhedron when it is created.
  static unsigned counter_;
};

std::ostream& operator<<(std::ostream& out, const Polyhedron& poly);


//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real Polyhedron::getDiameter() const
{
  return diameter_;
}

inline const Tetrahedron& Polyhedron::getTetra(SizeType i) const
{
  return tetrahedra_[i];
}

inline Tetrahedron& Polyhedron::getTetra(SizeType i)
{
  return tetrahedra_[i];
}

inline SizeType Polyhedron::getTetrahedraNo() const
{
  return tetrahedra_.size();
}

inline unsigned Polyhedron::getId() const
{
  return id_;
}

inline const Eigen::AlignedBox3d& Polyhedron::getBoundingBox() const
{
  return boundingBox_;
}

inline Polyhedron::ConstIterVertex Polyhedron::verticesCbegin() const
{
  return vertices_.cbegin();
}

inline Polyhedron::ConstIterVertex Polyhedron::verticesCend() const
{
  return vertices_.cend();
}

inline SizeType Polyhedron::getVerticesNo() const
{
  return vertices_.size();
}

inline void Polyhedron::resetCounter(unsigned counter)
{
  counter_ = counter;
}

} // namespace PolyDG

#endif // _POLYHEDRON_HPP_
