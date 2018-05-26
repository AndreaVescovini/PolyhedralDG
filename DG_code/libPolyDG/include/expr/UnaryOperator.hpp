/*!
    @file   UnaryOperator.hpp
    @author Andrea Vescovini
    @brief  Template class for unary operators
*/

#ifndef _UNARY_OPERATOR_HPP_
#define _UNARY_OPERATOR_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Template class for unary operators to be applied to an expression

    This template class is an expression and inherits from ExprWrapper<UnaryOperator<RO, OP>>.
    It defines a unary operator to be applied to an expression.

    @param RO The right operand.
    @param OP The operation.
*/

template <typename RO, typename OP>
class UnaryOperator : public ExprWrapper<UnaryOperator<RO, OP>>
{
public:
  //! Alias of return type of the unary operator
  using ReturnType = typename RO::ReturnType;

  //! Constructor
  explicit UnaryOperator(const RO& ro)
    : ro_{ro} {}

  //! Copy constructor
  UnaryOperator(const UnaryOperator&) = default;

  //! Move constructor
  UnaryOperator(UnaryOperator&&) = default;

  /*!
      @brief Call operator that evaluates the UnaryOperator inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return OP()(ro_(fe, i, t, p));
  }

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
  {
    return OP()(ro_(fe, i, j, t, p));
  }

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return OP()(ro_(fe, i, p));
  }

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
  {
    return OP()(ro_(fe, i, j, p));
  }

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceInt

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param sj Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
  {
    return OP()(ro_(fe, i, j, si, sj, p));
  }

  //! Destructor
  virtual ~UnaryOperator() = default;

private:
  //! Right operand
  const RO& ro_;
};

/*!
    @brief Negation

    This struct is a functor that defines the negation of a PolyDG::Real
    and of a @c Eigen::Vector3d.
*/

struct Negate
{
  //! Call operator for the negation of a PolyDG::Real
  Real operator()(Real ro) const
  {
    return -ro;
  }

  //! Call operator for the negation of a Eigen::Vector3d
  template <typename D>
  Eigen::Vector3d operator()(const Eigen::MatrixBase<D>& ro) const
  {
    return -ro;
  }
};

//! Overloading of the operator- for negating an expression
template <typename RO>
UnaryOperator<RO, Negate> operator-(const ExprWrapper<RO>& ro)
{
  return UnaryOperator<RO, Negate>(ro);
}

} // namespace PolyDG

#endif // _UNARY_OPERATOR_HPP_
