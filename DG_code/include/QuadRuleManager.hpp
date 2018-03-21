#ifndef _QUAD_RULE_MANAGER_HPP_
#define _QUAD_RULE_MANAGER_HPP_

#include "QuadRule.hpp"
#include <Eigen/Core>
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

// Const iterator over the vector of quadrature rules
  template <typename T>
  using CIter = typename std::vector<T>::const_iterator;

// Returns a quadrature rule over the standard tetrahedron of exactness doe.
// If it does not exist it returns the the quadrature rule with the maximum doe
// that is available.
  static const Rule3D& getTetraRule(unsigned doe);

// Returns a quadrature rule over the standard triangle of exactness doe.
// If it does not exist it returns the the quadrature rule with the maximum doe
// that is available.
  static const Rule2D& getTriaRule(unsigned doe);

// Returns the number of available quadrature rules over tetrahedra
  static unsigned getTetraRuleNo();

// Returns the number of available quadrature rules over triangles
  static unsigned getTriaRuleNo();

// Returns a const interator to the first quadrature rule over tetrahedra
  static CIter<Rule3D> tetraCbegin();

// Returns a const iterator to the rule after the end over tetrhedra
  static CIter<Rule3D> tetraCend();

// Returns a const interator to the first quadrature rule over triangles
  static CIter<Rule2D> triaCbegin();

// Returns a const iterator to the rule after the end over triangles
  static CIter<Rule2D> triaCend();

// Returns the matrix corresponding to the map from the standard 2d-simplex to the
// face faceNo of the standard 3d-simplex. The matrix is 3x3 so it must be applied
// to a 2d vector with homogeneous coordinates (i.e. a 3x1 vector with a 1 as
// last element)
  static const Eigen::Matrix3d& getFaceMap(unsigned faceNo);

  // static void addTetraRule(const QuadRule<Eigen::Vector3d>& rule);
  // static void addTriaRule(const QuadRule<Eigen::Vector2d>& rule);

private:
// Virtual constructor, copy-constructor and destructor, because this class
// contains only static methods and should be a simgleton
  QuadRuleManager() = default;
  QuadRuleManager(const QuadRuleManager& man) = default;
  virtual ~QuadRuleManager() = default;

// Vector containing the rules over the standard 3d-simplex
  static std::vector<Rule3D> tetraRules_;

// Vector containing the rules over the standard 2d-simplex
  static std::vector<Rule2D> triaRules_;

// Maps from the standard 2d-simplex to the four faces of the standard 3d-simplex
  static std::array<Eigen::Matrix3d, 4> faceMaps_;
};

}

#endif // _QUAD_RULE_MANAGER_HPP_
