/*!
    @file   PhiI.hpp
    @author Andrea Vescovini
    @brief  Class fot the expression for a test function
*/

#ifndef _PHI_I_
#define _PHI_I_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "PolyDG.hpp"

namespace PolyDG
{

/*!
    @brief Class fot the expression for a test function

    This class is an expression and inherits from ExprWrapper<PhiI>. It
    rapresents a test function \f$ \varphi_i \f$.
*/

class PhiI : public ExprWrapper<PhiI>
{
public:
  //! Alias for the return type of the call operator
  using ReturnType = Real;

  //! Constructor
  PhiI() = default;

  //! Copy constructor
  PhiI(const PhiI&) = default;

  //! Move constructor
  PhiI(PhiI&&) = default;

  /*!
      @brief Call operator that evaluates phi_i inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
  {
    return fe.getPhi(t, p, i);
  }

  /*!
      @brief Call operator that evaluates phi_i inside a FeElement

      The third argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const
  {
    return fe.getPhi(t, p, i);
  }

  /*!
      @brief Call operator that evaluates phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  ReturnType operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
  {
    return fe.getPhi(p, i);
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
    return fe.getPhi(p, i);
  }

  //! Destructor
  virtual ~PhiI() = default;
};

} // namespace PolyDG

#endif // _PHI_I_
