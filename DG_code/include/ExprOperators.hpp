#ifndef _EXPR_OPERATORS_HPP_
#define _EXPR_OPERATORS_HPP_

#include "geom.hpp"
#include <functional>
#include "ExprWrapper.hpp"
#include "FeElement.hpp"

namespace dgfem
{

// Class for binary operators
template <typename LO, typename RO, typename OP>
class BinaryOperator : public ExprWrapper<BinaryOperator<LO, RO, OP>>
{
public:
  BinaryOperator(const LO& lo, const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  virtual ~BinaryOperator() = default;

private:
  const LO& lo_;
  const RO& ro_;
};

// Specialization for a left operation by a scalar
template <typename RO, typename OP>
class BinaryOperator<geom::real, RO, OP> : public ExprWrapper<BinaryOperator<geom::real, RO, OP>>
{
public:
  BinaryOperator(geom::real lo, const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

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

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    geom::real ro_;
};

// Specialization for a left operation by an int
template <typename RO, typename OP>
class BinaryOperator<int, RO, OP> : public ExprWrapper<BinaryOperator<int, RO, OP>>
{
public:
  BinaryOperator(int lo, const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    int lo_;
    const RO& ro_;
};

// Specialization for a right operation by an int
template <typename LO, typename OP>
class BinaryOperator<LO, int, OP> : public ExprWrapper<BinaryOperator<LO, int, OP>>
{
public:
  BinaryOperator(const LO& lo, int ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  virtual ~BinaryOperator() = default;

  private:
    const LO& lo_;
    int ro_;
};

// Class for unary operators
template <typename RO, typename OP>
class UnaryOperator : public ExprWrapper<UnaryOperator<RO, OP>>
{
public:
  explicit UnaryOperator(const RO& ro);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;

  virtual ~UnaryOperator() = default;

private:
  const RO& ro_;
};

using Add = std::plus<geom::real>;
using Subtract = std::minus<geom::real>;
using Multiply = std::multiplies<geom::real>;
using Divide = std::divides<geom::real>;
using Negate = std::negate<geom::real>;

// Overloading of the operator + for the sum
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Add> operator+(const LO& lo, const RO& ro)
{
  return BinaryOperator<LO, RO, Add>(lo, ro);
}

// Overloading of the operator - for the subtraction
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Subtract> operator-(const LO& lo, const RO& ro)
{
  return BinaryOperator<LO, RO, Subtract>(lo, ro);
}

// Overloading of the operator * for the multuplication
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Multiply> operator*(const LO& lo, const RO& ro)
{
  return BinaryOperator<LO, RO, Multiply>(lo, ro);
}

// Overloading of the operator / for the division
template <typename LO, typename RO>
inline BinaryOperator<LO, RO, Divide> operator/(const LO& lo, const RO& ro)
{
  return BinaryOperator<LO, RO, Divide>(lo, ro);
}

// Overloading of the unary operator - for the negation
template <typename RO>
inline UnaryOperator<RO, Negate> operator-(const RO& ro)
{
  return UnaryOperator<RO, Negate>(ro);
}

//-------------------------------IMPLEMENTATION---------------------------------

template <typename LO, typename RO, typename OP>
BinaryOperator<LO, RO, OP>::BinaryOperator(const LO& lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename RO, typename OP>
geom::real BinaryOperator<LO, RO, OP>::operator()(const FeElement& fe,
                                                  unsigned i, unsigned j,
                                                  unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_(fe, i, j, t, q));
}

template <typename RO, typename OP>
BinaryOperator<geom::real, RO, OP>::BinaryOperator(geom::real lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename RO, typename OP>
geom::real BinaryOperator<geom::real, RO, OP>::operator()(const FeElement& fe,
                                                          unsigned i, unsigned j,
                                                          unsigned t, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, t, q));
}

template <typename LO, typename OP>
BinaryOperator<LO, geom::real, OP>::BinaryOperator(const LO& lo, geom::real ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename OP>
geom::real BinaryOperator<LO, geom::real, OP>::operator()(const FeElement& fe,
                                                          unsigned i, unsigned j,
                                                          unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_);
}

template <typename RO, typename OP>
BinaryOperator<int, RO, OP>::BinaryOperator(int lo, const RO& ro)
  : lo_{lo}, ro_{ro} {}

template <typename RO, typename OP>
geom::real BinaryOperator<int, RO, OP>::operator()(const FeElement& fe,
                                                   unsigned i, unsigned j,
                                                   unsigned t, unsigned q) const
{
  return OP()(lo_, ro_(fe, i, j, t, q));
}

template <typename LO, typename OP>
BinaryOperator<LO, int, OP>::BinaryOperator(const LO& lo, int ro)
  : lo_{lo}, ro_{ro} {}

template <typename LO, typename OP>
geom::real BinaryOperator<LO, int, OP>::operator()(const FeElement& fe,
                                                   unsigned i, unsigned j,
                                                   unsigned t, unsigned q) const
{
  return OP()(lo_(fe, i, j, t, q), ro_);
}

template <typename RO, typename OP>
UnaryOperator<RO, OP>::UnaryOperator(const RO& ro)
  : ro_{ro} {}

template <typename RO, typename OP>
geom::real UnaryOperator<RO, OP>::operator()(const FeElement& fe,
                                             unsigned i, unsigned j,
                                              unsigned t, unsigned q) const
{
  return OP()(ro_(fe, i, j, t, q));
}

}

#endif // _EXPR_OPERATORS_HPP_
