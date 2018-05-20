/*!
    @file   QuadRuleManager.hpp
    @author Andrea Vescovini
    @brief  Singleton class used to store and handle quadrature rules
*/

#ifndef _QUAD_RULE_MANAGER_HPP_
#define _QUAD_RULE_MANAGER_HPP_

#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <set>

namespace PolyDG
{

/*!
    @brief Singleton class used to store and handle quadrature rules

    This class is used to store, handle and provide qudrature rules, so it is
    implemented as a singleton in order to have one and only one object that
    communicates with the other classes to provide qudrarure rules. The
    implemenation exploits the @a Mayer's @a trick, so that the construtor is not
    available and you must rely in the method instance() to access the
    QuadRuleManager object.

    By default the following QuadRule3D are implemented:
    @arg 1  point,  degree of exactness 1, [Q].
    @arg 4  points, degree of exactness 2, [Q].
    @arg 5  points, degree of exactness 3, [Q].
    @arg 11 points, degree of exactness 4, with one negative weight, [K].
    @arg 15 points, degree of exactness 5, [K].
    @arg 24 points, degree of exactness 6, [K].
    @arg 31 points, degree of exactness 7, with one negative weight, [K].
    @arg 45 points, degree of exactness 8, with one negative weight, [K].

    By default the following QuadRule2D are implemented:
    @arg 1  point,  degree of exactness 1, [D].
    @arg 3  points, degree of exactness 2, [D].
    @arg 4  points, degree of exactness 3, with one negative weight, [D].
    @arg 6  points, degree of exactness 4, [D].
    @arg 7  points, degree of exactness 5, [D].
    @arg 12 points, degree of exactness 6, [D].
    @arg 13 points, degree of exactness 7, with one negative weight [D].
    @arg 16 points, degree of exactness 8, [D].
    @arg 19 points, degree of exactness 9, [D].
    @arg 25 points, degree of exactness 10, [D].

    [D] Dunavant D.: High degree ecient symmetrical Gaussian quadrature
        rules for the triangle. @a International @a Journal @a for @a Numerical
        @a Methods @a in @a Engineering, 21, 1129-1148 (1985).@n
    [K] Keast P.: Moderate-degree tetrahedral quadrature formulas. @a Computer
        @a Methods @a in @a Applied @a Mechanics @a and @a Engineering, 55,
        339-348 (1986).@n
    [Q] Quarteroni A.: @a Numerical @a Models @a for @a Differential @a Problems.
        MS&A, Springer-Verlag Italia, Milan (2014).

    When you ask for a quadraure rule with a certain degree of exactness, you
    get the first quadrature rule with a degree of exactness greater or equal to
    the one that you required. If not present you get the quadrature rule with
    the maximum degee of exactness that is available.@n
    You can add a new rule with setTriaRule and setTetraRule, but only one rule
    for each degree of exactness is allowed, so if you add a rule with a degree
    of exactness that is already present, the old one is erased.

    The class provides also four maps from the reference triangle in 2D to the
    four faces of the reference tetrahedron in 3D. With getFaceMap you get a
    @c Eigen::Matrix3d that should be applied to a vector in \f$\mathbb{R}^2\f$ with homogeneous
    coordinates (so with a 1 appended as third component).
*/

class QuadRuleManager
{
public:
  //! Alias for a bidirectional const iterator over quadrature rules
  template <typename T>
  using ConstIter = typename std::set<T, std::less<T>>::const_iterator;

  //! Get a reference to the singleton object
  static QuadRuleManager& instance();

  //! Deleted copy constructor
  QuadRuleManager(const QuadRule3D&) = delete;

  //! Deleted copy-assigment operator
  QuadRuleManager& operator=(const QuadRuleManager&) = delete;

  //! Deleted move-constructor
  QuadRuleManager(QuadRuleManager&&) = delete;

  //! Deleted move-assigment operator
  QuadRuleManager& operator=(QuadRuleManager&) = delete;

  /*!
      @brief Get a QuadRule3D

      This function returns a quadrature rule over the referece tetrahedron with
      exactness at least doe. If not implemented, it returns the the quadrature
      rule with the maximum degree of exactness that is available.
  */
  const QuadRule3D& getTetraRule(unsigned doe) const;

  /*!
      @brief Get a QuadRule3D

      This function returns a quadrature rule over the referece triangle with
      exactness at least doe. If not implemented, it returns the the quadrature
      rule with the maximum degree of exactness that is available.
  */
  const QuadRule2D& getTriaRule(unsigned doe) const;

  //! Get the number of available quadrature rules over tetrahedra
  inline SizeType getTetraRuleNo() const;

  //! Get the number of available quadrature rules over triangles
  inline SizeType getTriaRuleNo() const;

  //! Get a ConstIter pointing to the first QuadRule3D
  inline ConstIter<QuadRule3D> tetraCbegin() const;

  //! Get a ConstIter pointing to the @a past-the-end QuadRule3D
  inline ConstIter<QuadRule3D> tetraCend() const;

  //! Get a ConstIter pointing to the first QuadRule2D
  inline ConstIter<QuadRule2D> triaCbegin() const;

  //! Get a ConstIter pointing to the @a past-the-end QuadRule2D
  inline ConstIter<QuadRule2D> triaCend() const;

  /*!
      @brief Get a map from the reference triangle to the faces of the reference tetrahedron

      This function returns a matrix corresponding to the map from the 2d-simplex
      to the face faceNo of the 3d-simplex. The matrix is 3x3 so it must be
      applied to a vector in \f$\mathbb{R}^2\f$ with homogeneous coordinates
      (i.e. a 3x1 vector with a 1 as last element).

      @param faceNo Number from 0 to 3 corresponding to the face of the reference
                    tetrahedron to which the reference triangle is mapped. In the
                    reference tetrahedron face 0 is on the plane \f$ z = 0 \f$,
                    face 1 is one the plane \f$ y = 0 \f$, face 2 is on the plane
                    \f$ x = 0 \f$ and face 3 is on the plane \f$ x + y + z = 1 \f$.

  */
  inline const Eigen::Matrix3d& getFaceMap(unsigned faceNo) const;

  /*!
      @brief Set a QuadRule3D

      This function takes a const reference to a QuadRule3D and stores it in the
      QuadRuleManager. If there is already a rule with the same degree of
      exactness, that rule is substituted.
  */
  void setTetraRule(const QuadRule3D& rule);

  /*!
      @brief Set a QuadRule3D

      This function takes a rvalue reference to a QuadRule3D and moves it in the
      QuadRuleManager. If there is already a rule with the same degree of
      exactness, that rule is substituted.
  */
  void setTetraRule(QuadRule3D&& rule);

  /*!
      @brief Set a QuadRule2D

      This function takes a const reference to a QuadRule2D and stores it in the
      QuadRuleManager. If there is already a rule with the same degree of
      exactness, that rule is substituted.
  */
  void setTriaRule(const QuadRule2D& rule);

  /*!
      @brief Set a QuadRule2D

      This function takes a rvalue reference to a QuadRule3D and moves it in the
      QuadRuleManager. If there is already a rule with the same degree of
      exactness, that rule is substituted.
  */
  void setTriaRule(QuadRule2D&& rule);

  //! Destructor
  virtual ~QuadRuleManager() = default;

private:
  //! Constructor
  QuadRuleManager();

  //! Set containing the rules over the standard 3d-simplex
  std::set<QuadRule3D, std::less<QuadRule3D>> tetraRules_;

  //! Set containing the rules over the standard 2d-simplex
  std::set<QuadRule2D, std::less<QuadRule2D>> triaRules_;

  //! Maps from the standard 2d-simplex to the four faces of the standard 3d-simplex
  std::array<Eigen::Matrix3d, 4> faceMaps_;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline SizeType QuadRuleManager::getTetraRuleNo() const
{
  return tetraRules_.size();
}

inline SizeType QuadRuleManager::getTriaRuleNo() const
{
  return triaRules_.size();
}

inline QuadRuleManager::ConstIter<QuadRule3D> QuadRuleManager::tetraCbegin() const
{
  return tetraRules_.cbegin();
}

inline QuadRuleManager::ConstIter<QuadRule3D> QuadRuleManager::tetraCend() const
{
  return tetraRules_.cend();
}

inline QuadRuleManager::ConstIter<QuadRule2D>  QuadRuleManager::triaCbegin() const
{
  return triaRules_.cbegin();
}

inline QuadRuleManager::ConstIter<QuadRule2D>  QuadRuleManager::triaCend() const
{
  return triaRules_.cend();
}

inline const Eigen::Matrix3d& QuadRuleManager::getFaceMap(unsigned faceNo) const
{
  return faceMaps_[faceNo];
}

} // namespace PolyDG

#endif // _QUAD_RULE_MANAGER_HPP_
