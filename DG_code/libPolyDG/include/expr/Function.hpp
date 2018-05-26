/*!
    @file   Function.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for a function from R^3 to R
*/

#ifndef _FUNCTION_HPP_
#define _FUNCTION_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

#include <functional>

namespace PolyDG
{

/*!
    @brief Class for the expression for a function from R^3 to R

    This class is an expression and inherits from ExprWrapper<Function>. It
    rapresents a function \f$ f: \mathbb{R}^3 \rightarrow \mathbb{R} \f$.
*/

class Function : public ExprWrapper<Function>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Alias for a function taking a @c Eigen::Vector3d and returning a PolyDG::Real
  using funR3R1 = std::function<Real (const Eigen::Vector3d&)>;

  //! Constructor
  explicit Function(const funR3R1& fun)
    : fun_{fun} {}

  //! Copy constructor
  Function(const Function&) = default;

  //! Move constructor
  Function(Function&&) = default;

  /*!
      @brief Call operator that evaluates the function inside a FeElement

      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeElement& fe, unsigned /* i */, SizeType t, SizeType p) const
  {
    return fun_(fe.getQuadPoint(t, p));
  }

  /*!
      @brief Call operator that evaluates the function inside a FeElement

      The second and third arguments are not used.

      @param fe FeElement over which the evaluation has to be done.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeElement& fe, unsigned /* i */, unsigned /* j */, SizeType t, SizeType p) const
  {
    return fun_(fe.getQuadPoint(t, p));
  }

  /*!
      @brief Call operator that evaluates the function inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned /* i */, SizeType p) const
  {
    return fun_(fe.getQuadPoint(p));
  }

  /*!
      @brief Call operator that evaluates the function inside a FeFaceExt

      The second and third arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType p) const
  {
    return fun_(fe.getQuadPoint(p));
  }

  /*!
      @brief Call operator that evaluates the function inside a FeFaceInt

      The second, third, fourth and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType p) const
  {
    return fun_(fe.getQuadPoint(p));
  }

  //! Destructor
  virtual ~Function() = default;

private:
  //! The function
  funR3R1 fun_;
};

} // namespace PolyDG

#endif // _FUNCTION_HPP_
