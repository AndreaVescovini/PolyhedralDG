#ifndef _QUAD_RULE_HPP_
#define _QUAD_RULE_HPP_

#include <Eigen/Dense>
#include <vector>
#include <initializer_list>
#include "geom.hpp"

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

}

#endif // _QUAD_RULE_HPP_
