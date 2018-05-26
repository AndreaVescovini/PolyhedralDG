/*!
    @file   JumpPhiI.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the jump of a test function across a face
*/

#ifndef _JUMP_PHI_I_
#define _JUMP_PHI_I_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the jump of a test function across a face

    This class is an expression and inherits from ExprWrapper<JumpPhiI>. It
    rapresents the jump of a test function across a face, i.e.
    \f$ [\varphi_i] = \varphi_i^+ \mathbf{n}^+ + \varphi_i^- \mathbf{n}^- \f$ over
    internal faces and \f$ [\varphi_i] = \varphi_i \mathbf{n} \f$ over external faces.
*/

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  JumpPhiI() = default;

  //! Copy constructor
  JumpPhiI(const JumpPhiI&) = default;

  //! Move constructor
  JumpPhiI(JumpPhiI&&) = default;

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return fe.getPhi(p, i) * fe.getNormal();
  }

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceExt

      The third argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
  {
    return fe.getPhi(p, i) * fe.getNormal();
  }

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceInt

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned i, SideType si, SizeType p) const
  {
    return fe.getPhi(si, p, i) * fe.getNormal() * (si == Out ? 1 : -1);
  }

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceInt

      The third and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const
  {
    return fe.getPhi(si, p, i) * fe.getNormal() * (si == Out ? 1 : -1);
  }

  //! Destructor
  virtual ~JumpPhiI() = default;
};

} // namespace PolyDG

#endif // _JUMP_PHI_I_
