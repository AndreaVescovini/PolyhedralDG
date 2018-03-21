#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "FeElement.hpp"
#include "geom.hpp"

namespace dgfem
{

struct Stiff
{
  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                               unsigned t, unsigned q) const
  {
    return fe.getPhiDer(t, q, j).transpose() * fe.getPhiDer(t, q, i);
  }
};

}

#endif // _OPERATORS_HPP_
