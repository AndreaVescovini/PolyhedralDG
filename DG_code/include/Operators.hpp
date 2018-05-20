/*!
    @file   Operators.hpp
    @author Andrea Vescovini
    @brief  Classes that define expressions for the operators of a DG variational formulation
*/

#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "ExprWrapper.hpp"
#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "PolyDG.hpp"

#include <functional>

namespace PolyDG
{

/*!
    @brief Stiffness operator expression

    This class is an expression and inherits from ExprWrapper<Striff>. It
    rapresents the stiffness operator \f$ \nabla \varphi_j \cdot \nabla
    \varphi_i \f$, where \f$ \varphi_j \f$ is related to the solution and \f$ \varphi_i \f$
    is related to the test function.@n
    Using this operator is equivalent to using @c dot(uGrad,vGrad) where @c u is of type
    GradPhiJ and @c v is of type GradPhiI.
*/

class Stiff : public ExprWrapper<Stiff>
{
public:
  //! Constructor
  Stiff() = default;

  //! Copy constructor
  Stiff(const Stiff&) = default;

  //! Move constructor
  Stiff(Stiff&&) = default;

  /*!
      @brief Call operator that evaluates the stiffness term inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  //! Destructor
  virtual ~Stiff() = default;
};

/*!
    @brief Mass operator expression

    This class is an expression and inherits from ExprWrapper<Mass>. It
    rapresents the mass operator \f$ \varphi_j \varphi_i \f$, where \f$ \varphi_j \f$
    is related to the solution and \f$ \varphi_i \f$ is related to the test function.@n
    Using this operator is equivalent to using @c u*v where @c u is of type
    PhiJ and @c v is of type PhiI.
*/

class Mass : public ExprWrapper<Mass>
{
public:
  //! Constructor
  Mass() = default;

  //! Copy constructor
  Mass(const Mass&) = default;

  //! Move constructor
  Mass(Mass&&) = default;

  /*!
      @brief Call operator that evaluates the mass term inside a FeElement

      @param fe FeElement over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  //! Destructor
  virtual ~Mass() = default;
};

/*!
    @brief Expression for the gradient of a test function

    This class is an expression and inherits from ExprWrapper<GradPhiI>. It
    rapresents the gradient of a test function \f$ \nabla \varphi_i \f$.
*/

class GradPhiI : public ExprWrapper<GradPhiI>
{
public:
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
  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const;

  /*!
      @brief Call operator that evaluates the gradient of phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  //! Destructor
  virtual ~GradPhiI() = default;
};

/*!
    @brief Expression for the gradient of a basis function related to the solution

    This class is an expression and inherits from ExprWrapper<GradPhiJ>. It
    rapresents the gradient of a basis function related to the solution
    \f$ \nabla \varphi_j \f$.
*/

class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  //! Constructor
  GradPhiJ() = default;

  //! Copy constructor
  GradPhiJ(const GradPhiJ&) = default;

  //! Move constructor
  GradPhiJ(GradPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates the gradient of phi_j inside a FeElement

      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const;

  //! Destructor
  virtual ~GradPhiJ() = default;
};

/*!
    @brief Expression for a test function

    This class is an expression and inherits from ExprWrapper<PhiI>. It
    rapresents a test function \f$ \varphi_i \f$.
*/

class PhiI : public ExprWrapper<PhiI>
{
public:
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
  inline Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;

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
  inline Real operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const;

  /*!
      @brief Call operator that evaluates phi_i inside a FeFaceExt

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  //! Destructor
  virtual ~PhiI() = default;
};

/*!
    @brief Expression for a basis function related to the solution

    This class is an expression and inherits from ExprWrapper<PhiJ>. It
    rapresents a basis function related to the solution \f$ \varphi_j \f$.
*/

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  //! Constructor
  PhiJ() = default;

  //! Copy constructor
  PhiJ(const PhiJ&) = default;

  //! Move constructor
  PhiJ(PhiJ&&) = default;

  /*!
      @brief Call operator that evaluates phi_j inside a FeElement

      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const;

  //! Destructor
  virtual ~PhiJ() = default;
};

/*!
    @brief Expression for the jump of a test function across a face

    This class is an expression and inherits from ExprWrapper<JumpPhiI>. It
    rapresents the jump of a test function across a face, i.e.
    \f$ [\varphi_i] = \varphi_i^+ \mathbf{n}^+ + \varphi_i^- \mathbf{n}^- \f$ over
    internal faces and \f$ [\varphi_i] = \varphi_i \mathbf{n} \f$ over external faces.
*/

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  //! Constructor
  JumpPhiI() = default;

  //! Copy constructor
  JumpPhiI(const JumpPhiI&) = default;

  //! Move constructor
  JumpPhiI(JumpPhiI&&) = default;

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceInt

      The third and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const;

  /*!
      @brief Call operator that evaluates jump_phi_i inside a FeFaceExt

      The third argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const;

  //! Destructor
  virtual ~JumpPhiI() = default;
};

/*!
    @brief Expression for the jump of a basis function related to the solution across a face

    This class is an expression and inherits from ExprWrapper<JumpPhiJ>. It
    rapresents the jump of a basis function related to the solution across a face,
    i.e. \f$ [\varphi_j] = \varphi_j^+ \mathbf{n}^+ + \varphi_j^- \mathbf{n}^- \f$ over
    internal faces and \f$ [\varphi_j] = \varphi_j \mathbf{n} \f$ over external faces.
*/

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  //! Constructor
  JumpPhiJ() = default;

  //! Copy constructor
  JumpPhiJ(const JumpPhiJ&) = default;

  //! Move constructor
  JumpPhiJ(JumpPhiJ&&) = default;

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
  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates jump_phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~JumpPhiJ() = default;
};

/*!
    @brief Expression for the average of the gradient test function across a face

    This class is an expression and inherits from ExprWrapper<AverGradPhiI>. It
    rapresents the average of the gradient of a test function across a face, i.e.
    \f$ \{\!\!\{ \nabla \varphi_i \}\!\!\} = 0.5(\nabla \varphi_i^+ + \nabla \varphi_i^-) \f$ over
    internal faces and \f$ \{\!\!\{ \nabla \varphi_i \}\!\!\} = \nabla \varphi_i \f$
    over external faces.
*/

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  //! Constructor
  AverGradPhiI() = default;

  //! Copy constructor
  AverGradPhiI(const AverGradPhiI&) = default;

  //! Move constructor
  AverGradPhiI(AverGradPhiI&&) = default;

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceInt

      The third and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param si Side from which the evaluation of the test function has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const;

  /*!
      @brief Call operator that evaluates aver_grad_phi_i inside a FeFaceExt

      The third argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param i  Index related to the test function, it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const;

  //! Destructor
  virtual ~AverGradPhiI() = default;
};

/*!
    @brief Expression for the average of a basis function related to the solution across a face

    This class is an expression and inherits from ExprWrapper<AverGradPhiJ>. It
    rapresents the average of a basis function related to the solution across a face,
    i.e. \f$ \{\!\!\{ \nabla \varphi_j \}\!\!\} = 0.5(\nabla \varphi_j^+ + \nabla \varphi_j^-) \f$ over
    internal faces and \f$ \{\!\!\{ \nabla \varphi_j \}\!\!\} = \nabla \varphi_j \f$
    over external faces.
*/

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  //! Constructor
  AverGradPhiJ() = default;

  //! Copy constructor
  AverGradPhiJ(const AverGradPhiJ&) = default;

  //! Move constructor
  AverGradPhiJ(AverGradPhiJ&&) = default;

  /*!
      @brief Call operator that evaluates aver_grad_phi_j inside a FeFaceInt

      The second and fourth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param sj Side from which the evaluation of the basis function related to
                the solution has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const;

  /*!
      @brief Call operator that evaluates aver_grad_phi_j inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param j  Index related to the basis function related to the solution,
                it can be 0,...,fe.getDof() -1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const;

  //! Destructor
  virtual ~AverGradPhiJ() = default;
};

/*!
    @brief Expression for the penalty scaling over faces

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
  /*!
      @brief Constructor
      @param sigma Constant \f$ \sigma \f$ that multiplies the penalty parameter.
  */
  explicit PenaltyScaling(Real sigma = 1.0);

  //! Copy constructor
  PenaltyScaling(const PenaltyScaling&) = default;

  //! Move constructor
  PenaltyScaling(PenaltyScaling&&) = default;

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceExt

      The second and third arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  inline Real operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const;

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceInt

      Only the first argument is used.

      @param fe FeFaceInt over which the evaluation has to be done.
  */
  inline Real operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType /* p */) const;

  /*!
      @brief Call operator that evaluates the penalty parameter inside a FeFaceExt

      The second, third and fourth arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  inline Real operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const;

  //! Destructor
  virtual ~PenaltyScaling() = default;

private:
  //! The constant
  Real sigma_;
};

/*!
    @brief Expression for the unitary normal vector related to a FeFaceExt

    This class is an expression and inherits from ExprWrapper<Normal>. It
    rapresents unitary normal vector \f$ \mathbf{n} \f$ related to a FeFaceExt.
*/

class Normal : public ExprWrapper<Normal>
{
public:
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
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const;

  /*!
      @brief Call operator that evaluates the normal vector of a FeFaceExt

      The second, third and fourth arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
  */
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const;

  //! Destructor
  virtual ~Normal() = default;
};

/*!
    @brief Expression for a function from R^3 to R

    This class is an expression and inherits from ExprWrapper<Function>. It
    rapresents a function \f$ f: \mathbb{R}^3 \rightarrow \mathbb{R} \f$.
*/

class Function : public ExprWrapper<Function>
{
public:
  //! Alias for a function taking a @c Eigen::Vector3d and returning a PolyDG::Real
  using fun3real = std::function<Real (const Eigen::Vector3d&)>;

  //! Constructor
  explicit Function(const fun3real& fun);

  //! Copy constructor
  Function(const Function&) = default;

  //! Move constructor
  Function(Function&&) = default;

  /*!
      @brief Call operator that evaluates the function inside a FeElement

      The second argument is not used.

      @param fe FeElement over which the evaluation has to be done.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeElement& fe, unsigned /* i */, SizeType t, SizeType p) const;

  /*!
      @brief Call operator that evaluates the function inside a FeElement

      The second and third arguments are not used.

      @param fe FeElement over which the evaluation has to be done.
      @param t  Index related to the tetrahedron over which the evaluation has
                to be done, it can be 0,...,fe.getTetrahedraNo() - 1.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeElement& fe, unsigned /* i */, unsigned /* j */, SizeType t, SizeType p) const;

  /*!
      @brief Call operator that evaluates the function inside a FeFaceInt

      The second, third, fourth and fifth arguments are not used.

      @param fe FeFaceInt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType p) const;

  /*!
      @brief Call operator that evaluates the function inside a FeFaceExt

      The second argument is not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeFaceExt& fe, unsigned /* i */, SizeType p) const;

  /*!
      @brief Call operator that evaluates the function inside a FeFaceExt

      The second and third arguments are not used.

      @param fe FeFaceExt over which the evaluation has to be done.
      @param p  Index related to the quadrature point at which the evaluation
                has to be done, it can be 0,...,fe.getQuadPointsNo() - 1.
  */
  inline Real operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType p) const;

  //! Destructor
  virtual ~Function() = default;

private:
  //! The function
  fun3real fun_;
};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real Stiff::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return fe.getPhiDer(t, p, i).dot(fe.getPhiDer(t, p, j));
}

inline Real Mass::operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const
{
  return fe.getPhi(t, p, i) * fe.getPhi(t, p, j);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
{
  return fe.getPhiDer(t, p, i);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const
{
  return fe.getPhiDer(t, p, i);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return fe.getPhiDer(p, i);
}

inline const Eigen::Vector3d& GradPhiJ::operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const
{
  return fe.getPhiDer(t, p, j);
}

inline Real PhiI::operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const
{
  return fe.getPhi(t, p, i);
}

inline Real PhiI::operator()(const FeElement& fe, unsigned i, unsigned /* j */, SizeType t, SizeType p) const
{
  return fe.getPhi(t, p, i);
}

inline Real PhiI::operator()(const FeFaceExt& fe, unsigned i, SizeType p) const
{
  return fe.getPhi(p, i);
}

inline Real PhiJ::operator()(const FeElement& fe, unsigned /* i */, unsigned j, SizeType t, SizeType p) const
{
  return fe.getPhi(t, p, j);
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const
{
  return fe.getPhi(si, p, i) * fe.getNormal() * (si == Out ? 1 : -1);
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
{
  return fe.getPhi(p, i) * fe.getNormal();
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const
{
  return fe.getPhi(sj, p, j) * fe.getNormal() * (sj == Out ? 1 : -1);
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
{
  return fe.getPhi(p, j) * fe.getNormal();
}

inline Eigen::Vector3d AverGradPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, SizeType p) const
{
  return 0.5 * fe.getPhiDer(si, p, i);
}

inline const Eigen::Vector3d& AverGradPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, SizeType p) const
{
  return fe.getPhiDer(p, i);
}

inline Eigen::Vector3d AverGradPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, SizeType p) const
{
  return 0.5 * fe.getPhiDer(sj, p, j);
}

inline const Eigen::Vector3d& AverGradPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, SizeType p) const
{
  return fe.getPhiDer(p, j);
}

inline Real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline Real PenaltyScaling::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType /* p */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline Real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline const Eigen::Vector3d& Normal::operator()(const FeFaceExt& fe, unsigned /* i */, SizeType /* p */) const
{
  return fe.getNormal();
}

inline const Eigen::Vector3d& Normal::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType /* p */) const
{
  return fe.getNormal();
}

inline Real Function::operator()(const FeElement& fe, unsigned /* i */, SizeType t, SizeType p) const
{
  return fun_(fe.getQuadPoint(t, p));
}

inline Real Function::operator()(const FeElement& fe, unsigned /* i */, unsigned /* j */, SizeType t, SizeType p) const
{
  return fun_(fe.getQuadPoint(t, p));
}

inline Real Function::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, SizeType p) const
{
  return fun_(fe.getQuadPoint(p));
}

inline Real Function::operator()(const FeFaceExt& fe, unsigned /* i */, SizeType p) const
{
  return fun_(fe.getQuadPoint(p));
}

inline Real Function::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, SizeType p) const
{
  return fun_(fe.getQuadPoint(p));
}

} // namespace PolyDG

#endif // _OPERATORS_HPP_
