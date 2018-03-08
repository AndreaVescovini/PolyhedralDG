#ifndef _FE_SPACE_HPP_
#define _FE_SPACE_HPP_

#include "Mesh.hpp"
#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "QuadRule.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>

namespace dgfem {

class FeSpace
{
public:
  using theMesh = geom::Mesh;

  FeSpace(theMesh& Th, unsigned order);

  void setOrder(unsigned order);
  unsigned getOrder() const;
  unsigned getDofNo() const;

  virtual ~FeSpace() = default;

private:
  const theMesh& Th_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
  std::vector<std::array<unsigned, 3>> basisComposition_;
  std::vector<FeElement> feElements_;
  std::vector<FeFaceInt> feFacesInt_;
  std::vector<FeFaceExt> feFacesExt_;
  const QuadRule<Eigen::Vector3d>& tetraRule_;
  const QuadRule<Eigen::Vector2d>& triaRule_;

  void integerComposition();
  void initialize();


};

}

#endif // _FE_SPACE_HPP_
