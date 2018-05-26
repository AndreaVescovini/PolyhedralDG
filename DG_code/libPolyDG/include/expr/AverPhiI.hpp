/*!
    @file   AverPhiI.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the average of a test function across a face
*/

#ifndef _AVER_PHI_I_
#define _AVER_PHI_I_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the expression for the average of a test function across a face

    This class is an expression and inherits from ExprWrapper<AverPhiI>. It
    rapresents the average of a test function across a face,
    i.e. \f$ \{ \varphi_i \} = 0.5 \varphi_i^+ \varphi_i^- \f$ over internal
    faces and \f$ \{ \varphi_i \} = \varphi_i \f$ over external faces.
*/

class AverPhiI : public ExprWrapper<AverPhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  AverPhiI() = default;

  //! Copy constructor
  AverPhiI(const AverPhiI&) = default;

  //! Move constructor
  AverPhiI(AverPhiI&&) = default;

  /*!
      @brief Call operator that evaluates aver_phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return fe.getPhi(p, i);
  }

  /*!
      @brief Call operator that evaluates aver_phi_i inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  Real operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
  {
    return fe.getPhi(p, i);
  }

  /*!
      @brief Call operator that evaluates aver_phi_i inside a FeFaceInt

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
    return 0.5 * fe.getPhi(si, p, i);
  }

  /*!
      @brief Call operator that evaluates aver_phi_i inside a FeFaceInt

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
    return 0.5 * fe.getPhi(si, p, i);
  }

  //! Destructor
  virtual ~AverPhiI() = default;
};

} // namespace PolyDG

#endif // _AVER_PHI_I_
