#ifndef _FE_ELEMENT_HPP_
#define _FE_ELEMENT_HPP_

#include "Polyhedron.hpp"
#include "geom.hpp"
#include "QuadRule.hpp"
#include <Eigen/Dense>
#include <vector>
#include <array>

namespace dgfem
{

class FeElement
{
public:
  using theElem = geom::Polyhedron;

  FeElement(const theElem& elem, unsigned order, unsigned dofNo,
            const std::vector<std::array<unsigned, 3>>& basisComposition,
            const QuadRule<Eigen::Vector3d>& tetraRule);

  unsigned getOrder() const;
  unsigned getDofNo() const;

  geom::real getPhi(unsigned t, unsigned p, unsigned f) const;
  const Eigen::Vector3d& getPhiDer(unsigned t, unsigned p, unsigned f) const;

  virtual ~FeElement() = default;

private:
  const theElem& elem_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
  const QuadRule<Eigen::Vector3d>& tetraRule_;

  std::vector<geom::real> phi_;
  std::vector<Eigen::Vector3d> phiDer_;

  unsigned sub2ind(unsigned t, unsigned p, unsigned f) const;
};

}

#endif // _FE_ELEMENT_HPP_
