/*!
    @file   AverGradPhiI.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the average of the gradient of a test function across a face
*/

#ifndef _AVER_GRAD_PHI_I_
#define _AVER_GRAD_PHI_I_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the average of the gradient of a test function across a face

    This class is an expression and inherits from ExprWrapper<AverGradPhiI>. It
    rapresents the average of the gradient of a test function across a face, i.e.
    \f$ \{\!\!\{ \nabla \varphi_i \}\!\!\} = 0.5(\nabla \varphi_i^+ + \nabla \varphi_i^-) \f$ over
    internal faces and \f$ \{\!\!\{ \nabla \varphi_i \}\!\!\} = \nabla \varphi_i \f$
    over external faces.
*/

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  AverGradPhiI() = default;

  //! Copy constructor
  AverGradPhiI(const AverGradPhiI&) = default;

  //! Move constructor
  AverGradPhiI(AverGradPhiI&&) = default;

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return fe.getPhiDer(p, i);
  }

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceExt

      The third argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
  {
    return fe.getPhiDer(p, i);
  }

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceInt

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned i, SideType si, SizeType p) const
  {
    return 0.5 * fe.getPhiDer(si, p, i);
  }

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceInt

      The third and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const
  {
    return 0.5 * fe.getPhiDer(si, p, i);
  }

  //! Destructor
  virtual ~AverGradPhiI() = default;
};

} // namespace PolyDG

#endif // _AVER_GRAD_PHI_I_
