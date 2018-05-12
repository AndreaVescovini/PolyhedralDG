#ifndef _QUAD_RULE_MANAGER_HPP_
#define _QUAD_RULE_MANAGER_HPP_

#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <set>
// #include <vector>

namespace PolyDG
{

class QuadRuleManager
{
public:
  // Const iterator over the vector of quadrature rules.
  template <typename T>
  using ConstIter = typename std::set<T, std::less<T>>::const_iterator;

  // Function that returns a reference to the singleton object.
  static QuadRuleManager& instance();

  // Deleted copy-constructor.
  QuadRuleManager(const QuadRule3D&) = delete;

  // Deleted copy-assigment operator.
  QuadRuleManager& operator=(const QuadRuleManager&) = delete;

  // Deleted move-constructor.
  QuadRuleManager(QuadRuleManager&&) = delete;

  // Deleted move-assigment operator.
  QuadRuleManager& operator=(QuadRuleManager&) = delete;

  // Returns a quadrature rule over the standard tetrahedron of exactness at least
  // doe. If it does not exist it returns the the quadrature rule with the maximum doe
  // that is available.
  const QuadRule3D& getTetraRule(unsigned doe) const;

  // Returns a quadrature rule over the standard triangle of exactness at least doe.
  // If it does not exist it returns the the quadrature rule with the maximum doe
  // that is available.
  const QuadRule2D& getTriaRule(unsigned doe) const;

  // Returns the number of available quadrature rules over tetrahedra.
  inline SizeType getTetraRuleNo() const;

  // Returns the number of available quadrature rules over triangles.
  inline SizeType getTriaRuleNo() const;

  // Returns a const interator to the first quadrature rule over tetrahedra.
  inline ConstIter<QuadRule3D> tetraCbegin() const;

  // Returns a const iterator to the rule after the end over tetrhedra.
  inline ConstIter<QuadRule3D> tetraCend() const;

  // Returns a const interator to the first quadrature rule over triangles.
  inline ConstIter<QuadRule2D> triaCbegin() const;

  // Returns a const iterator to the rule after the end over triangles.
  inline ConstIter<QuadRule2D> triaCend() const;

  // Returns the matrix corresponding to the map from the standard 2d-simplex to
  // the face faceNo of the standard 3d-simplex. The matrix is 3x3 so it must be
  // applied to a 2d vector with homogeneous coordinates (i.e. a 3x1 vector with
  // a 1 as last element).
  inline const Eigen::Matrix3d& getFaceMap(unsigned faceNo) const;

  // Function that takes a const reference to a QuadRule3D and stores it in the
  // class. If there is already a rule with the same degree of exactness, that
  // rule is substituted.
  void setTetraRule(const QuadRule3D& rule);

  // Function that takes a rvalue reference to a QuadRule3D and stores it in the
  // class. If there is already a rule with the same degree of exactness, that
  // rule is substituted.
  void setTetraRule(QuadRule3D&& rule);

  // Function that takes a const reference to a QuadRule2D and stores it in the
  // class. If there is already a rule with the same degree of exactness, that
  // rule is substituted.
  void setTriaRule(const QuadRule2D& rule);

  // Function that takes a rvalue reference to a QuadRule2D and stores it in the
  // class. If there is already a rule with the same degree of exactness, that
  // rule is substituted.
  void setTriaRule(QuadRule2D&& rule);

  // Default virtual destructor.
  virtual ~QuadRuleManager() = default;

private:
  // Constructor.
  QuadRuleManager();

  // Set containing the rules over the standard 3d-simplex.
  std::set<QuadRule3D, std::less<QuadRule3D>> tetraRules_;

  // Set containing the rules over the standard 2d-simplex.
  std::set<QuadRule2D, std::less<QuadRule2D>> triaRules_;

  // Maps from the standard 2d-simplex to the four faces of the standard 3d-simplex.
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
