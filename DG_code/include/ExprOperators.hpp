#ifndef _EXPR_OPERATORS_HPP_
#define _EXPR_OPERATORS_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "Operators.hpp"
#include "PolyDG.hpp"

#include <functional>

namespace PolyDG
{

// Defining the scalar product
struct DotProduct
{
  Real operator()(const Eigen::Vector3d& lo, const Eigen::Vector3d& ro) const
  {
    return lo.dot(ro);
  }

  Real operator()(Real lo, Real ro) const
  {
    return lo * ro;
  }
};

// Class for binary operators
template <typename LO, typename RO, typename OP>
class BinaryOperator : public ExprWrapper<BinaryOperator<LO, RO, OP>>
{
public:
  BinaryOperator(const LO& lo, const RO& ro);

  BinaryOperator(const BinaryOperator&) = default;
  BinaryOperator(BinaryOperator&&) = default;

  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~BinaryOperator() = default;

private:
  const LO& lo_;
  const RO& ro_;
};

// Specialization for a left operation by a scalar
template <typename RO, typename OP>
class BinaryOperator<Real, RO, OP> : public ExprWrapper<BinaryOperator<Real, RO, OP>>
{
public:
  BinaryOperator(Real lo, const RO& ro);

  BinaryOperator(const BinaryOperator&) = default;
  BinaryOperator(BinaryOperator&&) = default;

  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~BinaryOperator() = default;

  private:
    Real lo_;
    const RO& ro_;
};

// Specialization for a right operation by a scalar
template <typename LO, typename OP>
class BinaryOperator<LO, Real, OP> : public ExprWrapper<BinaryOperator<LO, Real, OP>>
{
public:
  BinaryOperator(const LO& lo, Real ro);

  BinaryOperator(const BinaryOperator&) = default;
  BinaryOperator(BinaryOperator&&) = default;

  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    Real ro_;
};

// Specialization for the scalar product
template <typename LO, typename RO>
class BinaryOperator<LO, RO, DotProduct> : public ExprWrapper<BinaryOperator<LO, RO, DotProduct>>
{
public:
  BinaryOperator(const LO& lo, const RO& ro);

  BinaryOperator(const BinaryOperator&) = default;
  BinaryOperator(BinaryOperator&&) = default;

  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    const RO& ro_;
};

// Class for unary operators
template <typename RO, typename OP>
class UnaryOperator : public ExprWrapper<UnaryOperator<RO, OP>>
{
public:
  explicit UnaryOperator(const RO& ro);

  UnaryOperator(const UnaryOperator&) = default;
  UnaryOperator(UnaryOperator&&) = default;

  Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~UnaryOperator() = default;

private:
  const RO& ro_;
};

// I exploit the functors of the standard library
using Add = std::plus<Real>;
using Subtract = std::minus<Real>;
using Multiply = std::multiplies<Real>;
using Divide = std::divides<Real>;
using Negate = std::negate<Real>;

// Overloading of the operator + for the sum
template <typename LO, typename RO>
BinaryOperator<LO, RO, Add> operator+(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

// Overloading of the operator - for the subtraction
template <typename LO, typename RO>
BinaryOperator<LO, RO, Subtract> operator-(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

// Overloading of the operator * for the multiplication
template <typename LO, typename RO>
BinaryOperator<LO, RO, Multiply> operator*(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Multiply>(lo, ro);
}

template <typename LO>
BinaryOperator<LO, Real, Multiply> operator*(const ExprWrapper<LO>& lo, Real ro)
{
  return BinaryOperator<LO, Real, Multiply>(lo, ro);
}

template <typename RO>
BinaryOperator<Real, RO, Multiply> operator*(Real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<Real, RO, Multiply>(lo, ro);
}

// Overloading of the operator / for the division
template <typename LO, typename RO>
BinaryOperator<LO, RO, Divide> operator/(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Divide>(lo, ro);
}

template <typename LO>
BinaryOperator<LO, Real, Divide> operator/(const ExprWrapper<LO>& lo, Real ro)
{
  return BinaryOperator<LO, Real, Divide>(lo, ro);
}

template <typename RO>
BinaryOperator<Real, RO, Divide> operator/(Real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<Real, RO, Divide>(lo, ro);
}

// Overloading of the unary operator - for the negation
template <typename RO>
UnaryOperator<RO, Negate> operator-(const ExprWrapper<RO>& ro)
{
  return UnaryOperator<RO, Negate>(ro);
}

// Scalar Product
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

template <typename LO, typename RO>
BinaryOperator<LO, RO, DotProduct>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return DotProduct()(lo_(fe, i, j, t, p), ro_(fe, i, j, t, p));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const
{
  return DotProduct()(lo_(fe, i, j, si, sj, p), ro_(fe, i, j, si, sj, p));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return DotProduct()(lo_(fe, i, p), ro_(fe, i, p));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const
{
  return DotProduct()(lo_(fe, i, j, p), ro_(fe, i, j, p));
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

} // namespace PolyDG

#endif // _EXPR_OPERATORS_HPP_
