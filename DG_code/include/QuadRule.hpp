#ifndef _QUAD_RULE_HPP_
#define _QUAD_RULE_HPP_

#include "PolyDG.hpp"

#include <Eigen/Core>

#include <vector>
#include <initializer_list>
#include <numeric>
#include <cmath>
// #include <iostream>

namespace PolyDG
{

template <typename T>
class QuadRule
{
public:
  QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<Real>& weights);
  QuadRule(unsigned doe, std::initializer_list<T> points, std::initializer_list<Real> weights);

  unsigned getDoe() const;
  unsigned getPointsNo() const;
  const T& getPoint(unsigned i) const;
  Real getWeight(unsigned i) const;

// Functions that returns true if the sum of the weights is the volume of the
// simplex up to a tolerance tol and returns false otherwise.
  bool checkRule(Real tol = 1e-10) const;

  virtual ~QuadRule() = default;

private:
// Degree of exactness od the quadrature formula
  unsigned doe_;

// Vector of quadrature points
  std::vector<T> points_;

// Vector of quadrature weigths
  std::vector<Real> weights_;
};

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
unsigned QuadRule<T>::getPointsNo() const
{
  return weights_.size();
}

template <typename T>
const T& QuadRule<T>::getPoint(unsigned i) const
{
  return points_[i];
}

template <typename T>
Real QuadRule<T>::getWeight(unsigned i) const
{
  return weights_[i];
}

template <typename T>
bool QuadRule<T>::checkRule(Real tol) const
{
  Real sum = std::accumulate(weights_.cbegin(), weights_.cend(), static_cast<Real>(0.0));
  // std::cout << "Sum of weigths = " << sum << std::endl;

  // I use the Gamma function in order to compute the factorial
  Real volume = 1. / std::tgamma(points_[0].size() + 1);

  if(std::abs(sum - volume) < tol)
    return true;
  else
    return false;
}

} // namespace PolyDG

#endif // _QUAD_RULE_HPP_
