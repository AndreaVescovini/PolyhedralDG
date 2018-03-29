#ifndef _EXPR_OPERATORS_HPP_
#define _EXPR_OPERATORS_HPP_

#include "geom.hpp"
#include <functional>
#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "Operators.hpp"

namespace dgfem
{

// I exploit the functors of the standard library
using Add = std::plus<geom::real>;
using Subtract = std::minus<geom::real>;
using Multiply = std::multiplies<geom::real>;
using Divide = std::divides<geom::real>;
using Negate = std::negate<geom::real>;

// Defining the scalar product
struct DotProduct
{
  inline geom::real operator()(const Eigen::Vector3d& lo, const Eigen::Vector3d& ro) const
  {
    return lo.dot(ro);
  }

  inline geom::real operator()(geom::real lo, geom::real ro) const
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

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~BinaryOperator() = default;

private:
  const LO& lo_;
  const RO& ro_;
};

// // Specialization for a left operation by a scalar
template <typename RO, typename OP>
class BinaryOperator<geom::real, RO, OP> : public ExprWrapper<BinaryOperator<geom::real, RO, OP>>
{
public:
  BinaryOperator(geom::real lo, const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    geom::real lo_;
    const RO& ro_;
};

// Specialization for a right operation by a scalar
template <typename LO, typename OP>
class BinaryOperator<LO, geom::real, OP> : public ExprWrapper<BinaryOperator<LO, geom::real, OP>>
{
public:
  BinaryOperator(const LO& lo, geom::real ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    geom::real ro_;
};

// Specialization for the scalar product
template <typename LO, typename RO>
class BinaryOperator<LO, RO, DotProduct> : public ExprWrapper<BinaryOperator<LO, RO, DotProduct>>
{
public:
  BinaryOperator(const LO& lo, const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

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

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~UnaryOperator() = default;

private:
  const RO& ro_;
};

// Overloading of the operator + for the sum
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Add> operator+(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

template <typename LO>
inline BinaryOperator<LO, geom::real, Add> operator+(const ExprWrapper<LO>& lo, geom::real ro)
{
  return BinaryOperator<LO, geom::real, Add>(lo, ro);
}

template <typename RO>
inline BinaryOperator<geom::real, RO, Add> operator+(geom::real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<geom::real, RO, Add>(lo, ro);
}

// Overloading of the operator - for the subtraction
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Subtract> operator-(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

template <typename LO>
inline BinaryOperator<LO, geom::real, Subtract> operator-(const ExprWrapper<LO>& lo, geom::real ro)
{
  return BinaryOperator<LO, geom::real, Subtract>(lo, ro);
}

template <typename RO>
inline BinaryOperator<geom::real, RO, Subtract> operator-(geom::real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<geom::real, RO, Subtract>(lo, ro);
}

// Overloading of the operator * for the multiplication
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Multiply> operator*(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Multiply>(lo, ro);
}

template <typename LO>
inline BinaryOperator<LO, geom::real, Multiply> operator*(const ExprWrapper<LO>& lo, geom::real ro)
{
  return BinaryOperator<LO, geom::real, Multiply>(lo, ro);
}

template <typename RO>
inline BinaryOperator<geom::real, RO, Multiply> operator*(geom::real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<geom::real, RO, Multiply>(lo, ro);
}

// Overloading of the operator / for the division
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Divide> operator/(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, Divide>(lo, ro);
}

template <typename LO>
inline BinaryOperator<LO, geom::real, Divide> operator/(const ExprWrapper<LO>& lo, geom::real ro)
{
  return BinaryOperator<LO, geom::real, Divide>(lo, ro);
}

template <typename RO>
inline BinaryOperator<geom::real, RO, Divide> operator/(geom::real lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<geom::real, RO, Divide>(lo, ro);
}

// Overloading of the unary operator - for the negation
template <typename RO>
inline UnaryOperator<RO, Negate> operator-(const ExprWrapper<RO>& ro)
{
  return UnaryOperator<RO, Negate>(ro);
}

// Scalar Product
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, DotProduct> dot(const ExprWrapper<LO>& lo, const ExprWrapper<RO>& ro)
{
  return BinaryOperator<LO, RO, DotProduct>(lo, ro);
}


//-------------------------------IMPLEMENTATION---------------------------------

template <typename LO, typename RO, typename OP>
BinaryOperator<LO, RO, OP>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO, typename OP>
geom::real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_(fe, i, j, t, q));
}

template <typename LO, typename RO, typename OP>
geom::real BinaryOperator<LO, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return OP()(lo_(fe, i, j, side1, side2, q), ro_(fe, i, j, side1, side2, q));
}

template <typename LO, typename RO, typename OP>
geom::real BinaryOperator<LO, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_(fe, i, j, q), ro_(fe, i, j, q));
}

template <typename RO, typename OP>
BinaryOperator<geom::real, RO, OP>::BinaryOperator(geom::real lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename RO, typename OP>
geom::real BinaryOperator<geom::real, RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, t, q));
}

template <typename RO, typename OP>
geom::real BinaryOperator<geom::real, RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, side1, side2, q));
}

template <typename RO, typename OP>
geom::real BinaryOperator<geom::real, RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, q));
}

template <typename LO, typename OP>
BinaryOperator<LO, geom::real, OP>::BinaryOperator(const LO& lo, geom::real ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename OP>
geom::real BinaryOperator<LO, geom::real, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_);
}

template <typename LO, typename OP>
geom::real BinaryOperator<LO, geom::real, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return OP()(lo_(fe, i, j, side1, side2, q), ro_);
}

template <typename LO, typename OP>
geom::real BinaryOperator<LO, geom::real, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(lo_(fe, i, j, q), ro_);
}

template <typename LO, typename RO>
BinaryOperator<LO, RO, DotProduct>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO>
geom::real BinaryOperator<LO, RO, DotProduct>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return DotProduct()(lo_(fe, i, t, q), ro_(fe, j, t, q));
}

template <typename LO, typename RO>
geom::real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return DotProduct()(lo_(fe, i, side1, q), ro_(fe, j, side2, q));
}

template <typename LO, typename RO>
geom::real BinaryOperator<LO, RO, DotProduct>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return DotProduct()(lo_(fe, i, q), ro_(fe, j, q));
}

template <typename RO, typename OP>
UnaryOperator<RO, OP>::UnaryOperator(const RO& ro)
  : ro_{ro} {}

template <typename RO, typename OP>
geom::real UnaryOperator<RO, OP>::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return OP()(ro_(fe, i, j, t, q));
}

template <typename RO, typename OP>
geom::real UnaryOperator<RO, OP>::operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const
{
  return OP()(ro_(fe, i, j, side1, side2, q));
}

template <typename RO, typename OP>
geom::real UnaryOperator<RO, OP>::operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const
{
  return OP()(ro_(fe, i, j, q));
}

} // namespace dgfem

#endif // _EXPR_OPERATORS_HPP_
