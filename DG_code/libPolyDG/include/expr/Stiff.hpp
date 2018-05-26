/*!
    @file   Stiff.hpp
    @author Andrea Vescovini
    @brief  Class for the stiffness operator expression
*/

#ifndef _STIFF_HPP_
#define _STIFF_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the stiffness operator expression

    This class is an expression and inherits from ExprWrapper<Striff>. It
    rapresents the stiffness operator \f$ \nabla \varphi_j \cdot \nabla
    \varphi_i \f$, where \f$ \varphi_j \f$ is related to the solution and \f$ \varphi_i \f$
    is related to the test function.@n
    Using this operator is equivalent to using @c dot(uGrad,vGrad) where @c u is of type
    GradPhiJ and @c v is of type GradPhiI.
*/

class Stiff : public ExprWrapper<Stiff>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  Stiff() = default;

  //! Copy constructor
  Stiff(const Stiff&) = default;

  //! Move constructor
  Stiff(Stiff&&) = default;

  /*!
      @brief Call operator that evaluates the stiffness term inside a FeElement

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
    return fe.getPhiDer(t, p, i).dot(fe.getPhiDer(t, p, j));
  }

  /*!
      @brief Call operator that evaluates the stiffness term inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, unsigned  j , SizeType p) const
  {
    return fe.getPhiDer(p, i).dot(fe.getPhiDer(p, j));
  }

  //! Destructor
  virtual ~Stiff() = default;
};

} // namespace PolyDG

#endif // _STIFF_HPP_
