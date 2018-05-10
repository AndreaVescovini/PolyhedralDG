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

template <typename T>
class QuadRule
{
public:
  QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<Real>& weights);
  QuadRule(unsigned doe, std::initializer_list<T> points, std::initializer_list<Real> weights);

  // Default copy-constructor.
  QuadRule(const QuadRule&) = default;

  // Default copy-assigment operator.
  QuadRule& operator=(const QuadRule&) = default;

  // Default move cosntructor.
  QuadRule(QuadRule&&) = default;

  // Default move-assigment operator.
  QuadRule& operator=(QuadRule&&) = default;

  unsigned getDoe() const;
  SizeType getPointsNo() const;
  const T& getPoint(SizeType i) const;
  Real getWeight(SizeType i) const;

  // Functions that returns true if the sum of the weights is the volume of the
  // simplex up to a tolerance tol and returns false otherwise.
  bool checkRuleWeights(Real tol = 1e-10) const;

  // Default virtual destructor.
  virtual ~QuadRule() = default;

private:
  // Degree of exactness od the quadrature formula
  unsigned doe_;

  // Vector of quadrature points
  std::vector<T> points_;

  // Vector of quadrature weigths
  std::vector<Real> weights_;
};

using QuadRule3D = QuadRule<Eigen::Vector3d>;
using QuadRule2D = QuadRule<Eigen::Vector2d>;

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<Real>& weights)
  : doe_{doe}, points_{points}, weights_{weights} {}

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, std::initializer_list<T> points, std::initializer_list<Real> weights)
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

} // namespace PolyDG

#endif // _QUAD_RULE_HPP_
