#ifndef _EXPR_OPERATORS_HPP_
#define _EXPR_OPERATORS_HPP_

#include "PolyDG.hpp"
#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "Operators.hpp"

#include <functional>

namespace PolyDG
{

// I exploit the functors of the standard library
using Add = std::plus<Real>;
using Subtract = std::minus<Real>;
using Multiply = std::multiplies<Real>;
using Divide = std::divides<Real>;
using Negate = std::negate<Real>;

// using fun3 = std::function<Real (Eigen::Vector3d)>;

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

  Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~BinaryOperator() = default;

private:
  const LO& lo_;
  const RO& ro_;
};

// // Specialization for a left operation by a scalar
template <typename RO, typename OP>
class BinaryOperator<Real, RO, OP> : public ExprWrapper<BinaryOperator<Real, RO, OP>>
{
public:
  BinaryOperator(Real lo, const RO& ro);

  Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

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

  Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    Real ro_;
};

// // Specialization for a left operation by a function
// template <typename RO, typename OP>
// class BinaryOperator<fun3, RO, OP> : public ExprWrapper<BinaryOperator<fun3, RO, OP>>
// {
// public:
//   BinaryOperator(const fun3& lo, const RO& ro);
//
//   Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
//   Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
//   Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
//   Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
//   Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;
//
//   virtual ~BinaryOperator() = default;
//
//   private:
//     const fun3& lo_;
//     const RO& ro_;
// };
//
// // Specialization for a rigth operation by a function
// template <typename LO, typename OP>
// class BinaryOperator<LO, fun3, OP> : public ExprWrapper<BinaryOperator<LO, fun3, OP>>
// {
// public:
//   BinaryOperator(const LO& lo, const fun3& ro);
//
//   Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
//   Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
//   Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
//   Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
//   Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;
//
//   virtual ~BinaryOperator() = default;
//
//   private:
//     const LO& lo_;
//     const fun3& ro_;
// };

// Specialization for the scalar product
template <typename LO, typename RO>
class BinaryOperator<LO, RO, DotProduct> : public ExprWrapper<BinaryOperator<LO, RO, DotProduct>>
{
public:
  BinaryOperator(const LO& lo, const RO& ro);

  Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

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

  Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~UnaryOperator() = default;

private:
  const RO& ro_;
};

// Overloading of the operator + for the sum
template <typename LO, typename RO>
BinaryOperator<LO, RO, Add> operator+(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

// template <typename LO>
// BinaryOperator<LO, Real, Add> operator+(const ExprWrapper<LO>& lo, Real ro)
// {
//   return BinaryOperator<LO, Real, Add>(lo, ro);
// }
//
// template <typename RO>
// BinaryOperator<Real, RO, Add> operator+(Real lo, const ExprWrapper<RO>& ro)
// {
//   return BinaryOperator<Real, RO, Add>(lo, ro);
// }

// Overloading of the operator - for the subtraction
template <typename LO, typename RO>
BinaryOperator<LO, RO, Subtract> operator-(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

// template <typename LO>
// BinaryOperator<LO, Real, Subtract> operator-(const ExprWrapper<LO>& lo, Real ro)
// {
//   return BinaryOperator<LO, Real, Subtract>(lo, ro);
// }
//
// template <typename RO>
// BinaryOperator<Real, RO, Subtract> operator-(Real lo, const ExprWrapper<RO>& ro)
// {
//   return BinaryOperator<Real, RO, Subtract>(lo, ro);
// }

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
//
// template <typename LO>
// BinaryOperator<LO, fun3, Multiply> operator*(const ExprWrapper<LO>& lo, const fun3& ro)
// {
//   return BinaryOperator<LO, fun3, Multiply>(lo, ro);
// }
//
// template <typename RO>
// BinaryOperator<fun3, RO, Multiply> operator*(const fun3& lo, const ExprWrapper<RO>& ro)
// {
//   return BinaryOperator<fun3, RO, Multiply>(lo, ro);
// }

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
//
// template <typename LO>
// BinaryOperator<LO, fun3, Divide> operator/(const ExprWrapper<LO>& lo, const fun3& ro)
// {
//   return BinaryOperator<LO, fun3, Divide>(lo, ro);
// }
//
// template <typename RO>
// BinaryOperator<fun3, RO, Divide> operator/(const fun3& lo, const ExprWrapper<RO>& ro)
// {
//   return BinaryOperator<fun3, RO, Divide>(lo, ro);
// }

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


////////////////////////////////////////////////////////////////////////////////
//-------------------------------IMPLEMENTATION-------------------------------//
////////////////////////////////////////////////////////////////////////////////

template <typename LO, typename RO, typename OP>
BinaryOperator<LO, RO, OP>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, t, q), ro_(fe, i, t, q));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_(fe, i, j, t, q));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
{
  return OP()(lo_(fe, i, j, si, sj, q), ro_(fe, i, j, si, sj, q));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return OP()(lo_(fe, i, q), ro_(fe, i, q));
}

template <typename LO, typename RO, typename OP>
Real BinaryOperator<LO, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_(fe, i, j, q), ro_(fe, i, j, q));
}

template <typename RO, typename OP>
BinaryOperator<Real, RO, OP>::BinaryOperator(Real lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, t, q));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, t, q));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, si, sj, q));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, q));
}

template <typename RO, typename OP>
Real BinaryOperator<Real, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, q));
}

template <typename LO, typename OP>
BinaryOperator<LO, Real, OP>::BinaryOperator(const LO& lo, Real ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, t, q), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
{
  return OP()(lo_(fe, i, j, si, sj, q), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_(fe, i, j, q), ro_);
}

template <typename LO, typename OP>
Real BinaryOperator<LO, Real, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return OP()(lo_(fe, i, q), ro_);
}
//
// template <typename RO, typename OP>
// BinaryOperator<fun3, RO, OP>::BinaryOperator(const fun3& lo, const RO& ro)
//   : lo_{lo}, ro_{ro} {}
//
// template <typename RO, typename OP>
// Real BinaryOperator<fun3, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
// {
//   return OP()(lo_(fe.getQuadPoint(t, q)), ro_(fe, i, t, q));
// }
//
// template <typename RO, typename OP>
// Real BinaryOperator<fun3, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
// {
//   return OP()(lo_(fe.getQuadPoint(t, q)), ro_(fe, i, j, t, q));
// }
//
// template <typename RO, typename OP>
// Real BinaryOperator<fun3, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
// {
//   return OP()(lo_(fe.getQuadPoint(q)), ro_(fe, i, j, si, sj, q));
// }
//
// template <typename RO, typename OP>
// Real BinaryOperator<fun3, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
// {
//   return OP()(lo_(fe.getQuadPoint(q)), ro_(fe, i, q));
// }
//
// template <typename RO, typename OP>
// Real BinaryOperator<fun3, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
// {
//   return OP()(lo_(fe.getQuadPoint(q)), ro_(fe, i, j, q));
// }
//
// template <typename LO, typename OP>
// BinaryOperator<LO, fun3, OP>::BinaryOperator(const LO& lo, const fun3& ro)
//   : lo_{lo}, ro_{ro} {}
//
// template <typename LO, typename OP>
// Real BinaryOperator<LO, fun3, OP>::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
// {
//   return OP()(lo_(fe, i, t, q), ro_(fe.getQuadPoint(t, q)));
// }
//
// template <typename LO, typename OP>
// Real BinaryOperator<LO, fun3, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
// {
//   return OP()(lo_(fe, i, j, t, q), ro_(fe.getQuadPoint(t, q)));
// }
//
// template <typename LO, typename OP>
// Real BinaryOperator<LO, fun3, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
// {
//   return OP()(lo_(fe, i, j, si, sj, q), ro_(fe.getQuadPoint(q)));
// }
//
// template <typename LO, typename OP>
// Real BinaryOperator<LO, fun3, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
// {
//   return OP()(lo_(fe, i, q), ro_(fe.getQuadPoint(q)));
// }
//
// template <typename LO, typename OP>
// Real BinaryOperator<LO, fun3, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
// {
//   return OP()(lo_(fe, i, j, q), ro_(fe.getQuadPoint(q)));
// }

template <typename LO, typename RO>
BinaryOperator<LO, RO, DotProduct>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return DotProduct()(lo_(fe, i, j, t, q), ro_(fe, i, j, t, q));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
{
  return DotProduct()(lo_(fe, i, j, si, sj, q), ro_(fe, i, j, si, sj, q));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return DotProduct()(lo_(fe, i, q), ro_(fe, i, q));
}

template <typename LO, typename RO>
Real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return DotProduct()(lo_(fe, i, j, q), ro_(fe, i, j, q));
}

template <typename RO, typename OP>
UnaryOperator<RO, OP>::UnaryOperator(const RO& ro)
  : ro_{ro} {}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(ro_(fe, i, j, t, q));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const
{
  return OP()(ro_(fe, i, j, si, sj, q));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return OP()(ro_(fe, i, q));
}

template <typename RO, typename OP>
Real UnaryOperator<RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(ro_(fe, i, j, q));
}

} // namespace PolyDG

#endif // _EXPR_OPERATORS_HPP_
