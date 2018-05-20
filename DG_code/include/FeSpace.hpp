/*!
    @file   FeSpace.hpp
    @author Andrea Vescovini
    @brief  Class that defines finite elements spaces
*/

#ifndef _FE_SPACE_HPP_
#define _FE_SPACE_HPP_

#include "FeElement.hpp"
#include "FeFaceExt.hpp"
#include "FeFaceInt.hpp"
#include "Mesh.hpp"
#include "PolyDG.hpp"
#include "QuadRule.hpp"

#include <Eigen/Core>

#include <array>
#include <iostream>
#include <vector>

namespace PolyDG
{

/*!
    @brief Class that defines finite elements spaces

    This class defines a finite element space and stores FeElements, FeFacesExt
    and FeFacesInt.@n
    It rapresents the space:
    \f[
      \mathcal{D}_r(\mathcal{T}) = \{ v \in L^2(\Omega) : v|_\kappa \in
	    \mathbb{P}_r(\kappa) \quad \forall \kappa \in \mathcal{T}  \}
    \f]
    where \f$ \mathcal{T} \f$ is the Mesh, \f$ r \f$ is the degree of the space,
    \f$ \Omega \f$ is the domain, \f$ \kappa \f$ are the elements and
    \f$ \mathbb{P}_r(\kappa) \f$ is the standard space of polynomials of degree
    less then or equal to \f$ r \f$ over the element \f$ \kappa \f$.
    The basis functions are built in the following way:@n
    @arg Construct a Cartesian bounding box \f$ B_\kappa=I_1\times I_2 \times
         I_3 \quad \forall\kappa\in \mathcal{T} \f$, such that \f$ \bar{\kappa}
         \subseteq \bar{B_\kappa} \f$.
    @arg Then on every bounding box \f$ B_\kappa \f$ define a standard
         polynomial space \f$ \mathbb{P}_r(B_\kappa) \f$, spanned by a set of
         basis functions \f$ \{ \phi_{i,\kappa} \}_{i=1}^{dim(\mathbb{P}_r(B_\kappa))} \f$.
    @arg Use tensor-product (scaled) Legendre polynomials i.e., denoting by
         \f$ L_r(x) \f$ the Legendre polynomial of degree \f$ r \f$ defined on
         the interval \f$ [-1, 1] \f$, the corresponding scaled Legendre polynomial
         on the interval \f$ I_b = [x_1, x_2] \f$ may be defined by:
         \f[
           L_r^{[I_b]} (x) = \frac{1}{\sqrt{h_b}} L_r \bigg( \frac{x-m_b}{h_b} \bigg),
         \f]
         where \f$ h_b = (x_2-x_1)/2 \f$ and \f$ m_b = (x_1+x_2)/2 \f$.
    @arg The basis on the box \f$ B_\kappa \f$ is given by:
         \f[
           \phi_{\kappa,i}(\mathbf{x}) = L_{r_1}^{[1]}(x)L_{r_2}^{[2]}(y)L_{r_3}^{[3]}(z), \quad
           r_1+r_2+r_3 \leq r, \quad r_k \geq 0, \quad k = 1,2,3.
         \f]
    @arg Finally the polynomial basis over the polyhedral element \f$ \kappa \f$
         is defined simply restricting the support of \f$ \phi_{\kappa, i}(\mathbf{x}),
         \; i=1,\dots,dim(\mathbb{P}_r(B_\kappa)) \f$ to \f$ \kappa \f$, i.e.
         choosing  \f$ \phi_{\kappa, i}|_\kappa (\mathbf{x}), \; i=1,\dots,dim(\mathbb{P}_r(B_\kappa)) \f$.
*/

class FeSpace
{
public:
  //! Alias for a random access const iterator over FeElement, FeFaceExt or FeFaceInt
  template <typename T>
  using ConstIter = typename std::vector<T>::const_iterator;

  /*!
      @brief Constructor with degrees of exactness

      This constructor is used to specify the degree of exactness for the
      quadrature formulas.

      @param Th        The Mesh over which the space is built.
      @param degree    The degree of the space and its basis functions.
      @param doeQuad3D The required degree of exactness for the quadrature rule
                       over tetrahedra.
      @param doeQuad2D The required degree of exactness for the quadrature rule
                       over triangles.
  */
  FeSpace(Mesh& Th, unsigned degree, unsigned doeQuad3D, unsigned doeQuad2D);

  /*!
      @brief Constructor

      This constructor sets the quadrature rules in such a way that the integrals
      of the Poisson problem are computed exactly. This means a degree of
      exactness of 2 * (degree - 1) for the quadrature rule over tetrhedra
      beacuse there is only the stiffness integral and a degree of exactness of
      2 * degree for the quadrature rule over triangles.

      @param Th        The Mesh over which the space is built.
      @param degree    The degree of the space and its basis functions.
  */
  FeSpace(Mesh& Th, unsigned degree);

  //! Copy constructor
  FeSpace(const FeSpace&) = default;

  //! Move constructor
  FeSpace(FeSpace&&) = default;

  //! Get the degree of the space
  inline unsigned getDegree() const;

  /*!
      @brief Get the number of degrees of freedom

      This function returns the number of degrees of freedom in every finite
      element being in 3D it is (degree+1)*(degree+2)*(degree+3)/(3!).
  */
  inline unsigned getDof() const;

  /*!
      @brief Get a FeElement

      This functions returns the i-th FeElement.

      @param i The index of the FeElement required, it can be 0,..,getFeElementsNo() - 1.
  */
  inline const FeElement& getFeElement(SizeType i) const;

  /*!
      @brief Get a FeFaceExt

      This functions returns the i-th FeFaceExt.

      @param i The index of the FeFaceExt required, it can be 0,..,getFeFacesExtNo() - 1.
  */
  inline const FeFaceExt& getFeFaceExt(SizeType i) const;

  /*!
      @brief Get a FeFaceInt

      This functions returns the i-th FeFaceInt.

      @param i The index of the FeFaceInt required, it can be 0,..,getFeFacesIntNo() - 1.
  */
  inline const FeFaceInt& getFeFaceInt(SizeType i) const;

  //! Get the number of FeElement
  inline SizeType getFeElementsNo() const;

  //! Get the number of FeFaceExt
  inline SizeType getFeFacesExtNo() const;

  //! Get the number of FeFaceInt
  inline SizeType getFeFacesIntNo() const;

  //! Get the Mesh over which the FeSpace is built
  inline const Mesh& getMesh() const;

  /*!
      @brief Get the composition of the basis functions into monomials

      This functions returns a @c std::vector containing for each basis function
      an @c std::array with the powers at which the three monomials appear.@n
      For example if degree = 2 it contains [2, 0, 0], [1, 1, 0], [1, 0, 1],
      [0, 2, 0], [0, 1, 1], [0, 0, 2] that are the six possibility of compose a
      polynomial of degree 2 with three monomials.
  */
  inline const std::vector<std::array<unsigned, 3>>& getBasisComposition() const;

  //! Get a ConstIter pointing to the first FeElement
  inline ConstIter<FeElement> feElementsCbegin() const;

  //! Get a ConstIter pointing to the @a past-the-end FeElement
  inline ConstIter<FeElement> feElementsCend() const;

  //! Get a ConstIter pointing to the first FeFaceExt
  inline ConstIter<FeFaceExt> feFacesExtCbegin() const;

  //! Get a ConstIter pointing to the @a past-the-end FeFaceExt
  inline ConstIter<FeFaceExt> feFacesExtCend() const;

  //! Get a ConstIter pointing to the first FeFaceInt
  inline ConstIter<FeFaceInt> feFacesIntCbegin() const;

  //! Get a ConstIter pointing to the @a past-the-end FeFaceInt
  inline ConstIter<FeFaceInt> feFacesIntCend() const;

  //! Prints general information about the FeSpace
  void printInfo(std::ostream& out = std::cout) const;

  //! Prints the values of the basis functions over the elements
  void printElemBasis(std::ostream& out = std::cout) const;

  //! Prints the values of the gradients of the basis functions over the elements
  void printElemBasisDer(std::ostream& out = std::cout) const;

  //! Prints the values of the basis functions over the external and internal faces
  void printFaceBasis(std::ostream& out = std::cout) const;

  //! Prints the values of the gradients of the basis functions over the external and internal faces
  void printFaceBasisDer(std::ostream& out = std::cout) const;

  //! Destructor
  virtual ~FeSpace() = default;

private:
  //! Mesh over which the FeSpace is built
  const Mesh& Th_;

  //! Degree of polynomials
  const unsigned degree_;

  //! Number of degrees of freedom that in 3D is dof = (degree+1)*(degree+2)*(degree+3)/6
  const unsigned dof_;

  //! Possible degrees of the monomials that multiplied togheter give polynomials of degree less or equal to degree_
  std::vector<std::array<unsigned, 3>> basisComposition_;

  //! Vector of FeElement
  std::vector<FeElement> feElements_;

  //! Vector of FeFaceExt
  std::vector<FeFaceExt> feFacesExt_;

  //! Vector of FeFaceInt
  std::vector<FeFaceInt> feFacesInt_;

  //! Quadrature rule over tetrahedra
  const QuadRule3D& tetraRule_;

  //! Quadrature rule over triangles
  const QuadRule2D& triaRule_;

  //! Auxiliary function that computes basisComposition_.
  void integerComposition();

  //! Auxiliary function that fills feElements_, feFacesInt_ and feFacesExt_
  void initialize();

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeSpace::getDegree() const
{
  return degree_;
}

inline unsigned FeSpace::getDof() const
{
  return dof_;
}

inline const FeElement& FeSpace::getFeElement(SizeType i) const
{
  return feElements_[i];
}

inline const FeFaceExt& FeSpace::getFeFaceExt(SizeType i) const
{
  return feFacesExt_[i];
}

inline const FeFaceInt& FeSpace::getFeFaceInt(SizeType i) const
{
  return feFacesInt_[i];
}

inline SizeType FeSpace::getFeElementsNo() const
{
  return feElements_.size();
}

inline SizeType FeSpace::getFeFacesExtNo() const
{
  return feFacesExt_.size();
}

inline SizeType FeSpace::getFeFacesIntNo() const
{
  return feFacesInt_.size();
}

inline const Mesh& FeSpace::getMesh() const
{
  return Th_;
}

inline const std::vector<std::array<unsigned, 3>>& FeSpace::getBasisComposition() const
{
  return basisComposition_;
}

inline FeSpace::ConstIter<FeElement> FeSpace::feElementsCbegin() const
{
  return feElements_.cbegin();
}

inline FeSpace::ConstIter<FeElement> FeSpace::feElementsCend() const
{
  return feElements_.cend();
}

inline FeSpace::ConstIter<FeFaceExt> FeSpace::feFacesExtCbegin() const
{
  return feFacesExt_.cbegin();
}

inline FeSpace::ConstIter<FeFaceExt> FeSpace::feFacesExtCend() const
{
  return feFacesExt_.cend();
}

inline FeSpace::ConstIter<FeFaceInt> FeSpace::feFacesIntCbegin() const
{
  return feFacesInt_.cbegin();
}

inline FeSpace::ConstIter<FeFaceInt> FeSpace::feFacesIntCend() const
{
  return feFacesInt_.cend();
}

} // namespace PolyDG

#endif // _FE_SPACE_HPP_
