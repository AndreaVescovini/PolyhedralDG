/*!
    @file   ExprOperators.hpp
    @author Andrea Vescovini
    @brief  Template classes for binary and unary operators
*/

#ifndef _EXPR_OPERATORS_HPP_
#define _EXPR_OPERATORS_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "Operators.hpp"
#include "PolyDG.hpp"

#include <functional>

namespace PolyDG
{

/*!
    @brief Template class for binary operators between two expressions

    This template class is an expression and inherits from ExprWrapper<BinaryOperator<LO, RO, OP>>.
    It defines a binary operator between two expressions.

    @param LO The left operand.
    @param RO The right operand.
    @param OP The operation.
*/

template <typename LO, typename RO, typename OP>
class BinaryOperator : public ExprWrapper<BinaryOperator<LO, RO, OP>>
{
public:
  //! Constructor
  BinaryOperator(const LO& lo, const RO& ro);

  //! Copy constructor
  BinaryOperator(const BinaryOperator&) = default;

  //! Move constructor
  BinaryOperator(BinaryOperator&&) = default;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

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
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~BinaryOperator() = default;

private:
  //! The left operand
  const LO& lo_;

  //! The right operand
  const RO& ro_;
};

/*!
    @brief Specialization of BinaryOperator for a left operation with a PolyDG::Real

    This template class is an expression and inherits from ExprWrapper<BinaryOperator<PolyDG::Real, RO, OP>>.
    It is a specialization of BinaryOperator for a left operation with a PolyDG::Real.

*/

template <typename RO, typename OP>
class BinaryOperator<Real, RO, OP> : public ExprWrapper<BinaryOperator<Real, RO, OP>>
{
public:
  //! Constructor
  BinaryOperator(Real lo, const RO& ro);

  //! Copy constructor
  BinaryOperator(const BinaryOperator&) = default;

  //! Move constructor
  BinaryOperator(BinaryOperator&&) = default;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

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
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~BinaryOperator() = default;

  private:
    //! Left operand
    Real lo_;

    //! Right operand
    const RO& ro_;
};

/*!
    @brief Specialization of BinaryOperator for a right operation with a PolyDG::Real

    This template class is an expression and inherits from ExprWrapper<BinaryOperator<LO, PolyDG::Real, OP>>.
    It is a specialization of BinaryOperator for a right operation with a PolyDG::Real.

*/

template <typename LO, typename OP>
class BinaryOperator<LO, Real, OP> : public ExprWrapper<BinaryOperator<LO, Real, OP>>
{
public:
  //! Constructor
  BinaryOperator(const LO& lo, Real ro);

  //! Copy constructor
  BinaryOperator(const BinaryOperator&) = default;

  //! Move constructor
  BinaryOperator(BinaryOperator&&) = default;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

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
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~BinaryOperator() = default;

  private:
    //! Left operand
    const LO& lo_;

    //! Right operand
    Real ro_;
};

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
  //! Constructor
  explicit UnaryOperator(const RO& ro);

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
  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

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
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  /*!
      @brief Call operator that evaluates the BinaryOperator inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~UnaryOperator() = default;

private:
  //! Right operand
  const RO& ro_;
};

/*!
    @brief Scalar product

    This struct is a functor that defines the scalar product between two
    @c Eigen::Vector3d and by extension also the product between two PolyDG::Real.
*/
struct DotProduct
{
  //! Call operator that evaluates the dot product between two @c Eigen::Vector3d
  inline Real operator()(const Eigen::Vector3d& lo, const Eigen::Vector3d& ro) const;

  //! Call operator that evaluates the dot product between two PolyDG::Real
  inline Real operator()(Real lo, Real ro) const;
};

//! Alias for the addition functor for PolyDG::Real
using Add = std::plus<Real>;

//! Alias for the subtraction functor for PolyDG::Real
using Subtract = std::minus<Real>;

//! Alias for the multiplication functor for PolyDG::Real
using Multiply = std::multiplies<Real>;

//! Alias for the division functor for PolyDG::Real
using Divide = std::divides<Real>;

//! Alias for the negation functor for PolyDG::Real
using Negate = std::negate<Real>;

//! Overloading of the operator+ for summing expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Add> operator+(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

//! Overloading of the operator- for subtracting expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Subtract> operator-(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

//! Overloading of the operator* for multiplying expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Multiply> operator*(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Multiply>(lo, ro);
}

//! Overloading of the operator* for multiplying a PoltDG::Real and an expression
template <typename LO>
BinaryOperator<LO, Real, Multiply> operator*(const ExprWrapper<LO>& lo, Real ro)
{
  return BinaryOperator<LO, Real, Multiply>(lo, ro);
}

//! Overloading of the operator* for multiplying an expression and a PoltDG::Real
template <typename RO>
BinaryOperator<Real, RO, Multiply> operator*(Real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<Real, RO, Multiply>(lo, ro);
}

//! Overloading of the operator/ for dividing expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Divide> operator/(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Divide>(lo, ro);
}

//! Overloading of the operator/ for dividing an expression by a PolyDG::Real
template <typename LO>
BinaryOperator<LO, Real, Divide> operator/(const ExprWrapper<LO>& lo, Real ro)
{
  return BinaryOperator<LO, Real, Divide>(lo, ro);
}

//! Overloading of the operator/ fordividing a PoltDG::Real by an expression
template <typename RO>
BinaryOperator<Real, RO, Divide> operator/(Real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<Real, RO, Divide>(lo, ro);
}

//! Overloading of the operator- for the negating an expression
template <typename RO>
UnaryOperator<RO, Negate> operator-(const ExprWrapper<RO>& ro)
{
  return UnaryOperator<RO, Negate>(ro);
}

//! Operator dot(.,.) for the scalar product
template <typename LO, typename RO>
BinaryOperator<LO, RO, DotProduct> dot(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, DotProduct>(lo, ro);
}

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

template <typename LO, typename RO, typename OP>
BinaryOperator<LO, RO, OP>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
{
  return OP()(lo_(fe, i, t, p), ro_(fe, i, t, p));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return OP()(lo_(fe, i, j, t, p), ro_(fe, i, j, t, p));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
{
  return OP()(lo_(fe, i, j, si, sj, p), ro_(fe, i, j, si, sj, p));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return OP()(lo_(fe, i, p), ro_(fe, i, p));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
{
  return OP()(lo_(fe, i, j, p), ro_(fe, i, j, p));
}

template <typename RO, typename OP>
BinaryOperator<Real, RO, OP>::BinaryOperator(Real lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
{
  return OP()(lo_, ro_(fe, i, t, p));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return OP()(lo_, ro_(fe, i, j, t, p));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
{
  return OP()(lo_, ro_(fe, i, j, si, sj, p));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return OP()(lo_, ro_(fe, i, p));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
{
  return OP()(lo_, ro_(fe, i, j, p));
}

template <typename LO, typename OP>
BinaryOperator<LO, Real, OP>::BinaryOperator(const LO& lo, Real ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
{
  return OP()(lo_(fe, i, t, p), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return OP()(lo_(fe, i, j, t, p), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
{
  return OP()(lo_(fe, i, j, si, sj, p), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
{
  return OP()(lo_(fe, i, j, p), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return OP()(lo_(fe, i, p), ro_);
}

template <typename RO, typename OP>
UnaryOperator<RO, OP>::UnaryOperator(const RO& ro)
  : ro_{ro} {}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return OP()(ro_(fe, i, j, t, p));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
{
  return OP()(ro_(fe, i, j, si, sj, p));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return OP()(ro_(fe, i, p));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
{
  return OP()(ro_(fe, i, j, p));
}

inline Real DotProduct::operator()(const Eigen::Vector3d& lo, const Eigen::Vector3d& ro) const
{
  return lo.dot(ro);
}

inline Real DotProduct::operator()(Real lo, Real ro) const
{
  return lo * ro;
}

} // namespace PolyDG

#endif // _EXPR_OPERATORS_HPP_
