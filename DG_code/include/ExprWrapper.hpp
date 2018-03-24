#ifndef _EXPR_WRAPPER_HPP_
#define _EXPR_WRAPPER_HPP_

#include "geom.hpp"
#include "FeElement.hpp"

namespace dgfem
{

// Wrapper class that encapsulates every expression. It defines the cast
template <typename E>
class ExprWrapper
{
public:
  ExprWrapper() = default;

  void operator=(const ExprWrapper<E>&) = delete;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  operator const E&() const;
  operator E&();

  const E& asDerived() const;
  E& asDerived();

  virtual ~ExprWrapper() = default;
};

//-------------------------------IMPLEMENTATION---------------------------------

template <typename E>
geom::real ExprWrapper<E>::operator()(const FeElement& fe, unsigned i, unsigned j,
                                      unsigned t, unsigned q) const
{
  return asDerived().operator()(fe, i, j, t, q);
}

template <typename E>
ExprWrapper<E>::operator const E&() const
{
  return static_cast<const E&>(*this);
}

template <typename E>
ExprWrapper<E>::operator E&()
{
  return static_cast<E&>(*this);
}

template <typename E>
const E& ExprWrapper<E>::asDerived() const
{
  return static_cast<const E&>(*this);
}

template <typename E>
E& ExprWrapper<E>::asDerived()
{
  return static_cast<E&>(*this);
}

}

#endif // _EXPR_WRAPPER_HPP_
