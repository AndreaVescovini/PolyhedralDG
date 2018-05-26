/*!
    @file   JumpGradPhiJ.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the jump of the gradient of a basis function related to the solution across a face
*/

#ifndef _JUMP_GRAD_PHI_J_
#define _JUMP_GRAD_PHI_J_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the expression for the jump of the gradient of a basis function related to the solution across a face

    This class is an expression and inherits from ExprWrapper<JumpGradPhiJ>. It
    rapresents the jump of the gradient of a basis function related to the solution
    across a face, i.e. \f$ [\![ \nabla \varphi_j ]\!] = \nabla \varphi_j^+ \cdot \mathbf{n}^+ + \nabla \varphi_j^- \cdot \mathbf{n}^- \f$
    over internal faces and \f$ [\![ \nabla \varphi_j ]\!] = \nabla \varphi_j \cdot \mathbf{n} \f$ over external faces.
*/

class JumpGradPhiJ : public ExprWrapper<JumpGradPhiJ>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  JumpGradPhiJ() = default;

  //! Copy constructor
  JumpGradPhiJ(const JumpGradPhiJ&) = default;

  //! Move constructor
  JumpGradPhiJ(JumpGradPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates jump_grad_phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
  {
    return fe.getPhiDer(p, j).dot(fe.getNormal());
  }

  /*!
      @brief Call operator that evaluates jump_grad_phi_j inside a FeFaceInt

      The second and fourth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param sj Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const
  {
    return fe.getPhiDer(sj, p, j).dot(fe.getNormal()) * (sj == Out ? 1 : -1);
  }

  //! Destructor
  virtual ~JumpGradPhiJ() = default;
};

} // namespace PolyDG

#endif // _JUMP_GRAD_PHI_J_
