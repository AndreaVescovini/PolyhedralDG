#ifndef _EXPR_WRAPPER_HPP_
#define _EXPR_WRAPPER_HPP_

namespace dgfem
{

template <typename E>
class ExprWrapper {
public:
  ExprWrapper() = default;

  void operator=(const ExprWrapper<E>&) = delete;

  operator const E&() const;
  operator E&();

  virtual ~ExprWrapper() = default;
};

//-------------------------------IMPLEMENTATION---------------------------------

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

}

#endif // _EXPR_WRAPPER_HPP_
