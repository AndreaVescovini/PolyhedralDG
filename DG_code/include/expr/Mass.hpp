/*!
    @file   Mass.hpp
    @author Andrea Vescovini
    @brief  Class for the mass operator expression
*/

#ifndef _MASS_HPP_
#define _MASS_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the mass operator expression

    This class is an expression and inherits from ExprWrapper<Mass>. It
    rapresents the mass operator \f$ \varphi_j \varphi_i \f$, where \f$ \varphi_j \f$
    is related to the solution and \f$ \varphi_i \f$ is related to the test function.@n
    Using this operator is equivalent to using @c u*v where @c u is of type
    PhiJ and @c v is of type PhiI.
*/

class Mass : public ExprWrapper<Mass>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  Mass() = default;

  //! Copy constructor
  Mass(const Mass&) = default;

  //! Move constructor
  Mass(Mass&&) = default;

  /*!
      @brief Call operator that evaluates the mass term inside a FeElement

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
    return fe.getPhi(t, p, i) * fe.getPhi(t, p, j);
  }

  /*!
      @brief Call operator that evaluates the mass term inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, unsigned  j , SizeType p) const
  {
    return fe.getPhi(p, i) * fe.getPhi(p, j);
  }

  //! Destructor
  virtual ~Mass() = default;
};

} // namespace PolyDG

#endif // _MASS_HPP_
