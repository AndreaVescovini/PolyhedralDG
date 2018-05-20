/*!
    @file   QuadRule.hpp
    @author Andrea Vescovini
    @brief  Template class that defines a quadrature rule over a simplex
*/

#ifndef _QUAD_RULE_HPP_
#define _QUAD_RULE_HPP_

#include "PolyDG.hpp"

#include <Eigen/Core>

#include <cmath>
#include <initializer_list>
#include <numeric>
#include <vector>

namespace PolyDG
{

/*!
    @brief Template class that defines a quadrature rule over a simplex

    This template class defines a quadrature rule over a N-dimensional simplex
    that specify using as a template parameter T a vector rapresenting a point
    in dimension N. @n
    The two-dimensional simplex is the triangle (0, 0), (1, 0), (0, 1).@n
    The three-dimensional simplex is the tetrahedron (0, 0, 0), (1, 0, 0),
    (0, 1, 0), (0, 0, 1).@n
    Comparison functions based on the degree of exactness are provided.

    @param T Small vector that contains a quadrature point, for example
             @c Eigen::Vector2d or @c Eigen::Vector3d.
*/

template <typename T>
class QuadRule
{
public:
  /*!
      @brief Constructor
      @param doe     Degree of exactness of the rule.
      @param points  @c std::vector containing the quadrature points.
      @param weights @c std::vector containing the quadrature weights.
  */
  QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<Real>& weights);

  /*!
      @brief Constructor
      @param doe     Degree of exactness of the rule.
      @param points  @c std::initializer_list containing the quadrature points.
      @param weights @c std::initializer_list containing the quadrature weights.
  */
  QuadRule(unsigned doe, const std::initializer_list<T>& points, const std::initializer_list<Real>& weights);

  //! Copy constructor
  QuadRule(const QuadRule&) = default;

  //! Copy-assigment operator
  QuadRule& operator=(const QuadRule&) = default;

  //! Move cosntructor
  QuadRule(QuadRule&&) = default;

  //! Move-assigment operator
  QuadRule& operator=(QuadRule&&) = default;

  //! Get the degree of exactness
  unsigned getDoe() const;

  //! Get the number of quadrature points
  SizeType getPointsNo() const;

  /*!
      @brief Get a quadrature point

      This functions returns the i-th quadrature point.

      @param i The index of the quadrature point required, it can be 0,..,getPointsNo() - 1.
  */
  const T& getPoint(SizeType i) const;

  /*!
      @brief Get a quadrature weight

      This functions returns the i-th quadrature weight.

      @param i The index of the quadrature weight required, it can be 0,..,getPointsNo() - 1.
  */
  Real getWeight(SizeType i) const;

  /*!
      @brief Check the sum of the weights

      This function checks if the sum of the weights is the equal to volume of
      the simplex, as it has to be.

      @param  tol Tolerance up to which two PolyDG::Real ar considere the same.
      @return @c true if the check is successful, @c false if it is not.
  */
  bool checkRuleWeights(Real tol = 1e-10) const;

  //! Destructor.
  virtual ~QuadRule() = default;

  /*!
      @brief Binary operator that compares degrees of exactness

      This binary operator compares the degree of exactness of a rule with given
      one.

      @param rule The quadrature rule to be compared.
      @param doe  The value to be compared.
      @return @c true if the degree of exactness of rule is less then doe,
              @c false otherwise.
  */
  template <typename Q>
  friend bool compareDoe(const Q& rule, unsigned doe);

private:
  //! Degree of exactness od the quadrature rule
  unsigned doe_;

  //! Vector of quadrature points
  std::vector<T> points_;

  //! Vector of quadrature weigths
  std::vector<Real> weights_;
};

//! Alias for a QuadRule in three dimensions
using QuadRule3D = QuadRule<Eigen::Vector3d>;

//! Alias for a QuadRule in two dimensions
using QuadRule2D = QuadRule<Eigen::Vector2d>;

} // namespace PolyDG

namespace std
{

  //! Specialization of the functor @c std::less for a PolyDG::QuadRule based on the comparison of degrees of exactness
  template<typename T>
  struct less<PolyDG::QuadRule<T>>
  {
    //! Call operator
    bool operator()(const PolyDG::QuadRule<T>& lhs, const PolyDG::QuadRule<T>& rhs)
    {
      return lhs.getDoe() < rhs.getDoe();
    }
  };
} // namespace std

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

namespace PolyDG
{

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<Real>& weights)
  : doe_{doe}, points_{points}, weights_{weights} {}

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, const std::initializer_list<T>& points, const std::initializer_list<Real>& weights)
  : doe_{doe}, points_{points}, weights_{weights} {}

template <typename T>
unsigned QuadRule<T>::getDoe() const
{
  return doe_;
}

template <typename T>
SizeType QuadRule<T>::getPointsNo() const
{
  return weights_.size();
}

template <typename T>
const T& QuadRule<T>::getPoint(SizeType i) const
{
  return points_[i];
}

template <typename T>
Real QuadRule<T>::getWeight(SizeType i) const
{
  return weights_[i];
}

template <typename T>
bool QuadRule<T>::checkRuleWeights(Real tol) const
{
  const Real sum = std::accumulate(weights_.cbegin(), weights_.cend(), static_cast<Real>(0.0));

  // I use the Gamma function in order to compute the factorial.
  const Real volume = 1. / std::tgamma(points_[0].size() + 1);

  return (std::abs(sum - volume) < tol) ? true : false;
}

template <typename Q>
bool compareDoe(const Q& rule, unsigned doe)
{
  return rule.doe_ < doe;
}

} // namespace PolyDG

#endif // _QUAD_RULE_HPP_
