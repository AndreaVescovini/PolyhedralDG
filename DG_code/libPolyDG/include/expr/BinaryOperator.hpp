/*!
    @file   BinaryOperator.hpp
    @author Andrea Vescovini
    @brief  Template class for binary operators
*/

#ifndef _BINARY_OPERATOR_HPP_
#define _BINARY_OPERATOR_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

#include <utility>

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
  //! Alias of return type of the binary operator
  using ReturnType = decltype(OP()(std::declval<typename LO::ReturnType>(),
                                   std::declval<typename RO::ReturnType>()));

  //! Constructor
  BinaryOperator(const LO& lo, const RO& ro)
    : lo_{lo}, ro_{ro} {}

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
  ReturnType operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return OP()(lo_(fe, i, t, p), ro_(fe, i, t, p));
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
    return OP()(lo_(fe, i, j, t, p), ro_(fe, i, j, t, p));
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
    return OP()(lo_(fe, i, p), ro_(fe, i, p));
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
    return OP()(lo_(fe, i, j, p), ro_(fe, i, j, p));
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
    return OP()(lo_(fe, i, j, si, sj, p), ro_(fe, i, j, si, sj, p));
  }

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
  //! Alias of return type of the binary operator
  using ReturnType = decltype(OP()(std::declval<Real>(),
                                   std::declval<typename RO::ReturnType>()));

  //! Constructor
  BinaryOperator(Real lo, const RO& ro)
    : lo_{lo}, ro_{ro} {}

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
  ReturnType operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return OP()(lo_, ro_(fe, i, t, p));
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
    return OP()(lo_, ro_(fe, i, j, t, p));
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
    return OP()(lo_, ro_(fe, i, p));
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
    return OP()(lo_, ro_(fe, i, j, p));
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
    return OP()(lo_, ro_(fe, i, j, si, sj, p));
  }

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
  //! Alias of return type of the binary operator
  using ReturnType = decltype(OP()(std::declval<typename LO::ReturnType>(),
                                   std::declval<Real>()));

  //! Constructor
  BinaryOperator(const LO& lo, Real ro)
    : lo_{lo}, ro_{ro} {}

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
  ReturnType operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return OP()(lo_(fe, i, t, p), ro_);
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
    return OP()(lo_(fe, i, j, t, p), ro_);
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
    return OP()(lo_(fe, i, p), ro_);
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
    return OP()(lo_(fe, i, j, p), ro_);
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
    return OP()(lo_(fe, i, j, si, sj, p), ro_);
  }

  //! Destructor
  virtual ~BinaryOperator() = default;

  private:
    //! Left operand
    const LO& lo_;

    //! Right operand
    Real ro_;
};

/*!
    @brief Addition

    This struct is a functor that defines the addition between two PolyDG::Real
    and between two @c Eigen::Vector3d.
*/

struct Add
{
  //! Call operator that evaluates the addition between two PolyDG::Real
  Real operator()(Real lo, Real ro) const
  {
    return lo + ro;
  }

  //! Call operator that evaluates the addition between two @c Eigen::Vector3d
  template <typename D>
  Eigen::Vector3d operator()(const Eigen::MatrixBase<D>& lo, const Eigen::MatrixBase<D>& ro) const
  {
    return lo + ro;
  }
};

/*!
    @brief Subtraction

    This struct is a functor that defines the subtraction between two PolyDG::Real
    and between two @c Eigen::Vector3d.
*/

struct Subtract
{
  //! Call operator that evaluates the subtraction between two PolyDG::Real
  Real operator()(Real lo, Real ro) const
  {
    return lo - ro;
  }

  //! Call operator that evaluates the subtraction between two @c Eigen::Vector3d
  template <typename D>
  Eigen::Vector3d operator()(const Eigen::MatrixBase<D>& lo, const Eigen::MatrixBase<D>& ro) const
  {
    return lo - ro;
  }
};

/*!
    @brief Multiplication

    This struct is a functor that defines the multiplication between two PolyDG::Real
    and between a PolyDG::Real and a @c Eigen::Vector3d.
*/

struct Multiply
{
  //! Call operator that evaluates the multiplication between two PolyDG::Real
  Real operator()(Real lo, Real ro) const
  {
    return lo * ro;
  }

  //! Call operator that evaluates the multiplication between a @c Eigen::Vector3d and a PolyDG::Real
  template <typename D>
  Eigen::Vector3d operator()(const Eigen::MatrixBase<D>& lo, Real ro) const
  {
    return lo * ro;
  }

  //! Call operator that evaluates the multiplication between a PolyDG::Real and a @c Eigen::Vector3d
  template <typename D>
  Eigen::Vector3d operator()(Real lo, const Eigen::MatrixBase<D>& ro) const
  {
    return lo * ro;
  }
};

/*!
    @brief Division

    This struct is a functor that defines the division between two PolyDG::Real
    and between a @c Eigen::Vector3d and a PolyDG::Real.
*/

struct Divide
{
  //! Call operator that evaluates the division between two PolyDG::Real
  Real operator()(Real lo, Real ro) const
  {
    return lo / ro;
  }

  //! Call operator that evaluates the division between a @c Eigen::Vector3d and a PolyDG::Real
  template <typename D>
  Eigen::Vector3d operator()(const Eigen::MatrixBase<D>& lo, Real ro) const
  {
    return lo / ro;
  }
};

/*!
    @brief Scalar product

    This struct is a functor that defines the scalar product between two
    @c Eigen::Vector3d and by extension also the product between two PolyDG::Real.
*/
struct DotProduct
{
  //! Call operator that evaluates the dot product between two @c Eigen::Vector3d
  template <typename D>
  Real operator()(const Eigen::MatrixBase<D>& lo, const Eigen::MatrixBase<D>& ro) const
  {
    return lo.dot(ro);
  }

  //! Call operator that evaluates the dot product between two PolyDG::Real
  Real operator()(Real lo, Real ro) const
  {
    return lo * ro;
  }
};

//! Overloading of the operator+ for summing two expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Add> operator+(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

//! Overloading of the operator- for subtracting two expressions
template <typename LO, typename RO>
BinaryOperator<LO, RO, Subtract> operator-(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

//! Overloading of the operator* for multiplying two expressions
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

//! Overloading of the operator/ for dividing two expressions
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

//! Overloading of the operator/ for dividing a PoltDG::Real by an expression
template <typename RO>
BinaryOperator<Real, RO, Divide> operator/(Real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<Real, RO, Divide>(lo, ro);
}

//! Operator dot(.,.) for the scalar product
template <typename LO, typename RO>
BinaryOperator<LO, RO, DotProduct> dot(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, DotProduct>(lo, ro);
}

} // namespace PolyDG

#endif // _BINARY_OPERATOR_HPP_
