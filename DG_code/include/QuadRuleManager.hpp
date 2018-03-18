#ifndef _QUAD_RULE_MANAGER_HPP_
#define _QUAD_RULE_MANAGER_HPP_

#include "QuadRule.hpp"
#include <Eigen/Dense>
#include <vector>

namespace dgfem
{

class QuadRuleManager
{
public:
  using Rule3D = QuadRule<Eigen::Vector3d>;
  using Rule2D = QuadRule<Eigen::Vector2d>;

  // template <typename T>
  // using RuleContainer = std::vector<T>;

  template <typename T>
  using CIter = typename std::vector<T>::const_iterator;

  static const Rule3D& getTetraRule(unsigned doe);
  static const Rule2D& getTriaRule(unsigned doe);

  static unsigned getTetraRuleNo();
  static unsigned getTriaRuleNo();

  static CIter<Rule3D> tetraCbegin();
  static CIter<Rule3D> tetraCend();
  static CIter<Rule2D> triaCbegin();
  static CIter<Rule2D> triaCend();

  static const Eigen::Matrix3d& getFaceMap(unsigned faceNo);
  // static void addTetraRule(const QuadRule<Eigen::Vector3d>& rule);
  // static void addTriaRule(const QuadRule<Eigen::Vector2d>& rule);

private:
  QuadRuleManager() = default;
  QuadRuleManager(const QuadRuleManager& man) = default;
  virtual ~QuadRuleManager() = default;

  static std::vector<Rule3D> tetraRules_;
  static std::vector<Rule2D> triaRules_;

  static std::array<Eigen::Matrix3d, 4> faceMaps_;
};

}

#endif // _QUAD_RULE_MANAGER_HPP_
