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

class FeSpace
{
public:
  template <typename T>
  using ConstIter = typename std::vector<T>::const_iterator;

  // Constructor that takes a mesh Th, the degree of polynomials for the space and
  // the required degrees of exactness for quadrature formulas.
  FeSpace(Mesh& Th, unsigned degree, unsigned quad3DDegree, unsigned quad2DDegree);

  // Constructor that takes a mesh Th and the degree of polynomials for the space,
  // the quadrature formulas are chosen to fit with degree.
  FeSpace(Mesh& Th, unsigned degree);

  // Default copy-constructor.
  FeSpace(const FeSpace&) = default;

  // Default move-constructor.
  FeSpace(FeSpace&&) = default;

  inline unsigned getdegree() const;
  inline unsigned getDof() const;

  inline const FeElement& getFeElement(SizeType i) const;
  inline const FeFaceExt& getFeFaceExt(SizeType i) const;
  inline const FeFaceInt& getFeFaceInt(SizeType i) const;

  inline SizeType getFeElementsNo() const;
  inline SizeType getFeFacesExtNo() const;
  inline SizeType getFeFacesIntNo() const;

  inline const Mesh& getMesh() const;
  inline const std::vector<std::array<unsigned, 3>>& getBasisComposition() const;

  inline ConstIter<FeElement> feElementsCbegin() const;
  inline ConstIter<FeElement> feElementsCend() const;
  inline ConstIter<FeFaceExt> feFacesExtCbegin() const;
  inline ConstIter<FeFaceExt> feFacesExtCend() const;
  inline ConstIter<FeFaceInt> feFacesIntCbegin() const;
  inline ConstIter<FeFaceInt> feFacesIntCend() const;

  void printInfo(std::ostream& out = std::cout) const;
  void printElemBasis(std::ostream& out = std::cout) const;
  void printElemBasisDer(std::ostream& out = std::cout) const;
  void printFaceBasis(std::ostream& out = std::cout) const;
  void printFaceBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeSpace() = default;

private:
  // Reference to the mesh over which the FeSpace is built.
  const Mesh& Th_;

  // degree of polynomials.
  const unsigned degree_;

  // Number of degrees of freedom that in 3D is dof = (degree+1)*(degree+2)*(degree+3)/6
  const unsigned dof_;

  std::vector<std::array<unsigned, 3>> basisComposition_;
  std::vector<FeElement> feElements_;
  std::vector<FeFaceExt> feFacesExt_;
  std::vector<FeFaceInt> feFacesInt_;
  const QuadRule3D& tetraRule_;
  const QuadRule2D& triaRule_;

  // Auxiliary function that computes basisComposition_.
  void integerComposition();

  // Auxiliary function that fills feElements_, feFacesInt_ and feFacesExt_.
  void initialize();

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeSpace::getdegree() const
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
