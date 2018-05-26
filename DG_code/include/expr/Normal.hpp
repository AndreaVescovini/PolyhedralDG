/*!
    @file   Normal.hpp
    @author Andrea Vescovini
    @brief  Class for the expression for the unitary normal vector related to a FeFaceExt
*/

#ifndef _NORMAL_HPP_
#define _NORMAL_HPP_

#include "ExprWrapper.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

#include <Eigen/Core>

namespace PolyDG
{

/*!
    @brief Class for the expression for the unitary normal vector related to a FeFaceExt

    This class is an expression and inherits from ExprWrapper<Normal>. It
    rapresents unitary normal vector \f$ \mathbf{n} \f$ related to a face.
*/

class Normal : public ExprWrapper<Normal>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Eigen::Vector3d;

  //! Constructor
  Normal() = default;

  //! Copy constructor
  Normal(const Normal&) = default;

  //! Move constructor
  Normal(Normal&&) = default;

  /*!
      @brief Call operator that evaluates the normal vector of a FeFaceExt

      The second and third arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const
  {
    return fe.getNormal();
  }

  /*!
      @brief Call operator that evaluates the normal vector of a FeFaceExt

      The second, third and fourth arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const
  {
    return fe.getNormal();
  }

  //! Destructor
  virtual ~Normal() = default;
};

} // namespace PolyDG

#endif // _NORMAL_HPP_
