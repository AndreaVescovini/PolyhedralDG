#ifndef _EXPR_HPP_
#define _EXPR_HPP_

#include"FeElement.hpp"

namespace dgfem
{

template <typename T>
class Expr {
public:
  Expr() = default;
  explicit Expr(const T& form);

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                               unsigned t, unsigned q) const;

  virtual ~Expr() = default;

private:
  T formulation_;
};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename T>
Expr<T>::Expr(const T& form)
  : formulation_{form} {}

template <typename T>
geom::real Expr<T>::operator()(const FeElement& fe, unsigned i, unsigned j,
                               unsigned t, unsigned q) const
{
  return formulation_(fe, i, j, t, q);
}

}

#endif // _EXPR_HPP_
