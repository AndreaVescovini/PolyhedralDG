/*!
    @file   AverGradPhiJ.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the average of the gradient of a basis function related to the solution across a face
*/

#ifndef _AVER_GRAD_PHI_J_
#define _AVER_GRAD_PHI_J_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the average of the gradient of a basis function related to the solution across a face

    This class is an expression and inherits from ExprWrapper<AverGradPhiJ>. It
    rapresents the average of the gradient of a basis function related to the
    solution across a face, i.e. \f$ \{\!\!\{ \nabla \varphi_j \}\!\!\} = 0.5(\nabla \varphi_j^+ + \nabla \varphi_j^-) \f$
    over internal faces and \f$ \{\!\!\{ \nabla \varphi_j \}\!\!\} = \nabla \varphi_j \f$
    over external faces.
*/

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  AverGradPhiJ() = default;

  //! Copy constructor
  AverGradPhiJ(const AverGradPhiJ&) = default;

  //! Move constructor
  AverGradPhiJ(AverGradPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates aver_grad_phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
  {
    return fe.getPhiDer(p, j);
  }

  /*!
      @brief Call operator that evaluates aver_grad_phi_j inside a FeFaceInt

      The second and fourth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param sj Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const
  {
    return 0.5 * fe.getPhiDer(sj, p, j);
  }

  //! Destructor
  virtual ~AverGradPhiJ() = default;
};

} // namespace PolyDG

#endif // _AVER_GRAD_PHI_J_
