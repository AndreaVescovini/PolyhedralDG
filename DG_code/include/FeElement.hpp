#ifndef _FE_ELEMENT_HPP_
#define _FE_ELEMENT_HPP_

#include "Polyhedron.hpp"
#include "geom.hpp"
#include <Eigen/Dense>
#include <vector>

namespace dgfem
{

template <typename T>
using multivector = std::vector<std::vector<std::vector<T>>>;

class FeElement
{
public:
  using theElem = geom::Polyhedron;

  explicit FeElement(const theElem& elem);

  virtual ~FeElement() = default;
private:
  const theElem& elem_;
  // Jacobian qui anzichè in tetra?
  multivector<geom::real> phi;
  multivector<Eigen::Vector3d> dphi;
};

}

#endif // _FE_ELEMENT_HPP_
