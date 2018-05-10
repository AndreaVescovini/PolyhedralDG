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

  // Default copy-constructor.
  ExprWrapper(const ExprWrapper&) = default;

  // Deleted copy-assignment operator.
  ExprWrapper& operator=(const ExprWrapper&) = delete;

  // Default move-constructor.
  ExprWrapper(ExprWrapper&&) = default;

  // Deleted move-assignament operator.
  ExprWrapper& operator=(ExprWrapper&&) = delete;

  // Const cast operators to derived.
  operator const E&() const;

  // Cast operator to derived.
  operator E&();

  // Explicit method for the const cast to derived.
  const E& asDerived() const;

  // Explicit method fot the cast to derived.
  E& asDerived();

  // Default virtual destructor.
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
