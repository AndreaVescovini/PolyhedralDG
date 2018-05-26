/*!
    @file   PenaltyScaling.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the penalty scaling over faces
*/

#ifndef _PENALTY_SCALING_HPP_
#define _PENALTY_SCALING_HPP_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class for the expression for the penalty scaling over faces

    This class is an expression and inherits from ExprWrapper<PenaltyScaling>.
    It rapresents the penalty scaling over faces, i.e.
    \f$ \sigma \max\limits_{\kappa \in \{\kappa^+, \kappa^-\}} \big\{ \frac{r^2}{h_\kappa}\big\} \f$
    over internal faces and \f$ \sigma\frac{r^2}{h_\kappa} \f$ over external faces,
    where \f$ r \f$ is the degree of the FeSpace, \f$ h_{\kappa} \f$ is the
    diameter of the element \f$ \kappa \f$ and \f$ \kappa^+, \kappa^- \f$ are
    the two elements @a "In" and @a "Out" sharing the internal face.
*/
class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  /*!
      @brief Constructor
      @param sigma Constant \f$ \sigma \f$ that multiplies the penalty parameter.
  */
  explicit PenaltyScaling(Real sigma = 1.0)
    : sigma_{sigma} {}

  //! Copy constructor
  PenaltyScaling(const PenaltyScaling&) = default;

  //! Move constructor
  PenaltyScaling(PenaltyScaling&&) = default;

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceExt

      The second and third arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  Real operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const
  {
    return sigma_ * fe.getPenaltyParam();
  }

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceExt

      The second, third and fourth arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  Real operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const
  {
    return sigma_ * fe.getPenaltyParam();
  }

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceInt

      Only the first argument is used.

      @param fe FeFaceInt over which the evaluation has to be done.
  */
  Real operator()(const FeFaceInt& fe, unsigned /* i */, SideType /* si */, SizeType /* p */) const
  {
    return sigma_ * fe.getPenaltyParam();
  }

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceInt

      Only the first argument is used.

      @param fe FeFaceInt over which the evaluation has to be done.
  */
  Real operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType /* p */) const
  {
    return sigma_ * fe.getPenaltyParam();
  }

  //! Destructor
  virtual ~PenaltyScaling() = default;

private:
  //! The constant
  Real sigma_;
};

} // namespace PolyDG

#endif // _PENALTY_SCALING_HPP_
