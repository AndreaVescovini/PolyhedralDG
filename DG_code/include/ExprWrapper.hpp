#ifndef _EXPR_WRAPPER_HPP_
#define _EXPR_WRAPPER_HPP_

#include "geom.hpp"
// #include "FeElement.hpp"
// #include "FeFaceInt.hpp"
// #include "FeFaceExt.hpp"

namespace dgfem
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

// // Call operator for expressions that involve volume integrals
//   geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
//                         unsigned t, unsigned q) const;
//
// // Call operator for expressions that involve volume integrals, returning a vector
//   Eigen::Vector3d operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
//
// // Call operator for expressions that involve integrals over internal faces
//   geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j,
//                         int side1, int side2, unsigned q) const;
//
// // Call operator for expressions that involve integrals over external faces
//   geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;


// Cast operators
  operator const E&() const;
  operator E&();

// Explicit methods for the casting
  // const E& asDerived() const;
  // E& asDerived();

  virtual ~ExprWrapper() = default;
};

//-------------------------------IMPLEMENTATION---------------------------------

// template <typename E>
// geom::real ExprWrapper<E>::operator()(const FeElement& fe, unsigned i, unsigned j,
//                                       unsigned t, unsigned q) const
// {
//   return asDerived().operator()(fe, i, j, t, q);
// }
//
// template <typename E>
// Eigen::Vector3d ExprWrapper<E>::operator()(const FeElement& fe, unsigned i,
//                                            unsigned t, unsigned q) const
// {
//   return asDerived().operator()(fe, i, t, q);
// }
//
// template <typename E>
// geom::real ExprWrapper<E>::operator()(const FeFaceInt& fe, unsigned i, unsigned j,
//                                       int side1, int side2, unsigned q) const
// {
//   return asDerived().operator()(fe, i, j, side1, side2, q);
// }
//
// template <typename E>
// geom::real ExprWrapper<E>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
// {
//   return asDerived().operator()(fe, i, j, q);
// }

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

// template <typename E>
// const E& ExprWrapper<E>::asDerived() const
// {
//   return static_cast<const E&>(*this);
// }
//
// template <typename E>
// E& ExprWrapper<E>::asDerived()
// {
//   return static_cast<E&>(*this);
// }

}

#endif // _EXPR_WRAPPER_HPP_
