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
  static const QuadRule<Eigen::Vector3d>& getTetraRule(unsigned doe);
  static const QuadRule<Eigen::Vector2d>& getTriaRule(unsigned doe);

  // static void addTetraRule(const QuadRule<Eigen::Vector3d>& rule);
  // static void addTriaRule(const QuadRule<Eigen::Vector2d>& rule);

private:
  QuadRuleManager() = default;
  QuadRuleManager(const QuadRuleManager& man) = default;
  virtual ~QuadRuleManager() = default;

  static std::vector<QuadRule<Eigen::Vector3d>> tetraRules_;
  static std::vector<QuadRule<Eigen::Vector2d>> triaRules_;

  static std::vector<QuadRule<Eigen::Vector3d>> initTetraRules();
  static std::vector<QuadRule<Eigen::Vector2d>> initTriaRules();
};

}

#endif // _QUAD_RULE_MANAGER_HPP_
