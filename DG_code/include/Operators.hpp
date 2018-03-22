#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "FeElement.hpp"
#include "geom.hpp"
#include "ExprWrapper.hpp"

namespace dgfem
{

class Stiff : public ExprWrapper<Stiff>
{
public:
  Stiff() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;
  virtual ~Stiff() = default;

// private:
  // geom::real viscosity;

};

}

#endif // _OPERATORS_HPP_
