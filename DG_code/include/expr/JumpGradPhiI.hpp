/*!
    @file   JumpGradPhiI.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the jump of the gradient of a test function across a face
*/

#ifndef _JUMP_GRAD_PHI_I_
#define _JUMP_GRAD_PHI_I_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the expression for the jump of the gradient of a test function across a face

    This class is an expression and inherits from ExprWrapper<JumpGradPhiI>. It
    rapresents the jump of the gradient of a basis function related to the solution
    across a face, i.e. \f$ [\![ \nabla \varphi_i ]\!] = \nabla \varphi_i^+ \cdot \mathbf{n}^+ + \nabla \varphi_i^- \cdot \mathbf{n}^- \f$
    over internal faces and \f$ [\![ \nabla \varphi_i ]\!] = \nabla \varphi_i \cdot \mathbf{n} \f$ over external faces.
*/

class JumpGradPhiI : public ExprWrapper<JumpGradPhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  JumpGradPhiI() = default;

  //! Copy constructor
  JumpGradPhiI(const JumpGradPhiI&) = default;

  //! Move constructor
  JumpGradPhiI(JumpGradPhiI&&) = default;

  /*!
      @brief Call operator that evaluates jump_grad_phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return fe.getPhiDer(p, i).dot(fe.getNormal());
  }

  /*!
      @brief Call operator that evaluates jump_grad_phi_i inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
  {
    return fe.getPhiDer(p, i).dot(fe.getNormal());
  }

  /*!
      @brief Call operator that evaluates jump_grad_phi_i inside a FeFaceInt

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceInt& fe, unsigned i, SideType si, SizeType p) const
  {
    return fe.getPhiDer(si, p, i).dot(fe.getNormal()) * (si == Out ? 1 : -1);
  }

  /*!
      @brief Call operator that evaluates jump_grad_phi_i inside a FeFaceInt

      The second and fourth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const
  {
    return fe.getPhiDer(si, p, i).dot(fe.getNormal()) * (si == Out ? 1 : -1);
  }

  //! Destructor
  virtual ~JumpGradPhiI() = default;
};

} // namespace PolyDG

#endif // _JUMP_GRAD_PHI_I_
