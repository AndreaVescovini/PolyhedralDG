/*!
    @file   GradPhiI.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the gradient of a test function
*/

#ifndef _GRAD_PHI_I_
#define _GRAD_PHI_I_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the gradient of a test function

    This class is an expression and inherits from ExprWrapper<GradPhiI>. It
    rapresents the gradient of a test function \f$ \nabla \varphi_i \f$.
*/

class GradPhiI : public ExprWrapper<GradPhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  GradPhiI() = default;

  //! Copy constructor
  GradPhiI(const GradPhiI&) = default;

  //! Move constructor
  GradPhiI(GradPhiI&&) = default;

  /*!
      @brief Call operator that evaluates the gradient of phi_i inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return fe.getPhiDer(t, p, i);
  }

  /*!
      @brief Call operator that evaluates the gradient of phi_i inside a FeElement

      The third argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const
  {
    return fe.getPhiDer(t, p, i);
  }

  /*!
      @brief Call operator that evaluates the gradient of phi_i inside a FeFaceExt

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
      @brief Call operator that evaluates the gradient of phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  const ReturnType& operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
  {
    return fe.getPhiDer(p, i);
  }

  //! Destructor
  virtual ~GradPhiI() = default;
};

} // namespace PolyDG

#endif // _GRAD_PHI_I_
