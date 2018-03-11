#ifndef _QUAD_RULE_HPP_
#define _QUAD_RULE_HPP_

#include <Eigen/Dense>
#include <vector>
#include <initializer_list>
#include <numeric>
#include <cmath>
#include "geom.hpp"
// #include <iostream>

namespace dgfem
{

template <typename T>
class QuadRule
{
public:
  QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<geom::real>& weights);
  QuadRule(unsigned doe, std::initializer_list<T> points, std::initializer_list<geom::real> weights);

  unsigned getDoe() const;
  unsigned getPointsNo() const;
  const T& getPoint(unsigned i) const;
  geom::real getWeight(unsigned i) const;

  bool checkRule(geom::real tol = 1e-10) const;

  virtual ~QuadRule() = default;

private:
  unsigned doe_; // degree of exactness
  std::vector<T> points_;
  std::vector<geom::real> weights_;
};

//------------IMPLEMENTATION----------------------------------------------------

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, const std::vector<T>& points, const std::vector<geom::real>& weights)
  : doe_{doe}, points_{points}, weights_{weights} {}

template <typename T>
QuadRule<T>::QuadRule(unsigned doe, std::initializer_list<T> points, std::initializer_list<geom::real> weights)
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
geom::real QuadRule<T>::getWeight(unsigned i) const
{
  return weights_[i];
}

template <typename T>
bool QuadRule<T>::checkRule(geom::real tol) const
{
  geom::real sum = std::accumulate(weights_.cbegin(), weights_.cend(), static_cast<geom::real>(0.0));
  // std::cout << "Sum of weigths = " << sum << std::endl;
  geom::real volume = 1. / std::tgamma(points_[0].size() + 1);

  if(std::abs(sum - volume) < tol)
    return true;
  else
    return false;
}

}

#endif // _QUAD_RULE_HPP_
