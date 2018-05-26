/*!
    @file   JumpPhiJ.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the jump of a basis function related to the solution across a face
*/

#ifndef _JUMP_PHI_J_
#define _JUMP_PHI_J_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the jump of a basis function related to the solution across a face

    This class is an expression and inherits from ExprWrapper<JumpPhiJ>. It
    rapresents the jump of a basis function related to the solution across a face,
    i.e. \f$ [\varphi_j] = \varphi_j^+ \mathbf{n}^+ + \varphi_j^- \mathbf{n}^- \f$ over
    internal faces and \f$ [\varphi_j] = \varphi_j \mathbf{n} \f$ over external faces.
*/

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  JumpPhiJ() = default;

  //! Copy constructor
  JumpPhiJ(const JumpPhiJ&) = default;

  //! Move constructor
  JumpPhiJ(JumpPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates jump_phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
  {
    return fe.getPhi(p, j) * fe.getNormal();
  }

  /*!
      @brief Call operator that evaluates jump_phi_j inside a FeFaceInt

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
    return fe.getPhi(sj, p, j) * fe.getNormal() * (sj == Out ? 1 : -1);
  }

  //! Destructor
  virtual ~JumpPhiJ() = default;
};

} // namespace PolyDG

#endif // _JUMP_PHI_J_
