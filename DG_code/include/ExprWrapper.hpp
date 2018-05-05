#ifndef _EXPR_WRAPPER_HPP_
#define _EXPR_WRAPPER_HPP_

#include "PolyDG.hpp"

namespace PolyDG
{

// Wrapper class that encapsulates every expression. It defines the cast operator
// to the template parameter class, i.e. to all those class that will inherit from
// ExprWrapper.
template <typename E>
class ExprWrapper
{
public:
  ExprWrapper() = default;

  void operator=(const ExprWrapper<E>&) = delete;

// Cast operators
  operator const E&() const;
  operator E&();

// Explicit methods for the casting
  const E& asDerived() const;
  E& asDerived();

  virtual ~ExprWrapper() = default;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename E>
ExprWrapper<E>::operator const E&() const
{
  return *static_cast<const E*>(this);
}

template <typename E>
ExprWrapper<E>::operator E&()
{
  return *static_cast<E*>(this);
}

template <typename E>
const E& ExprWrapper<E>::asDerived() const
{
  return *static_cast<const E*>(this);
}

template <typename E>
E& ExprWrapper<E>::asDerived()
{
  return *static_cast<E*>(this);
}

} // namespace PolyDG

#endif // _EXPR_WRAPPER_HPP_
