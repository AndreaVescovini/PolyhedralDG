/*!
    @file   PhiJ.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for a basis function related to the solution
*/

#ifndef _PHI_J_
#define _PHI_J_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the expression for a basis function related to the solution

    This class is an expression and inherits from ExprWrapper<PhiJ>. It
    rapresents a basis function related to the solution \f$ \varphi_j \f$.
*/

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  PhiJ() = default;

  //! Copy constructor
  PhiJ(const PhiJ&) = default;

  //! Move constructor
  PhiJ(PhiJ&&) = default;

  /*!
      @brief Call operator that evaluates phi_j inside a FeElement

      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const
  {
    return fe.getPhi(t, p, j);
  }

  /*!
      @brief Call operator that evaluates _phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
  {
    return fe.getPhi(p, j);
  }

  //! Destructor
  virtual ~PhiJ() = default;
};

} // namespace PolyDG

#endif // _PHI_J_
