/*!
    @file   GradPhiJ.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the gradient of a basis function related to the solution
*/

#ifndef _GRAD_PHI_J_
#define _GRAD_PHI_J_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the gradient of a basis function related to the solution

    This class is an expression and inherits from ExprWrapper<GradPhiJ>. It
    rapresents the gradient of a basis function related to the solution
    \f$ \nabla \varphi_j \f$.
*/

class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  GradPhiJ() = default;

  //! Copy constructor
  GradPhiJ(const GradPhiJ&) = default;

  //! Move constructor
  GradPhiJ(GradPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates the gradient of phi_j inside a FeElement


      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const
  {
    return fe.getPhiDer(t, p, j);
  }

  /*!
      @brief Call operator that evaluates the gradient of phi_j inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
  {
    return fe.getPhiDer(p, j);
  }

  //! Destructor
  virtual ~GradPhiJ() = default;
};

} // namespace PolyDG

#endif // _GRAD_PHI_J_
