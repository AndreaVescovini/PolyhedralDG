#ifndef _FE_SPACE_HPP_
#define _FE_SPACE_HPP_

#include "Mesh.hpp"
#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "QuadRuleManager.hpp"

#include <Eigen/Core>

#include <vector>
#include <array>
#include <iostream>

namespace PolyDG
{

class FeSpace
{
public:
  using TheMesh = PolyDG::Mesh;

// Constructor that takes a mesh Th, the order of polynomials order and the required
// degrees of exactness for quadrature formulas
  FeSpace(TheMesh& Th, unsigned order, unsigned quad3DDegree, unsigned quad2DDegree);

// Constructor that takes a mesh Th and the order of polynomials order, the quadrature
// formulas are chosen to fit with order.
  FeSpace(TheMesh& Th, unsigned order);

  void setOrder(unsigned order);
  inline unsigned getOrder() const;
  inline unsigned getDofNo() const;

  inline const FeElement& getFeElement(unsigned i) const;
  inline const FeFaceInt& getFeFaceInt(unsigned i) const;
  inline const FeFaceExt& getFeFaceExt(unsigned i) const;

  inline unsigned getFeElementsNo() const;
  inline unsigned getFeFacesIntNo() const;
  inline unsigned getFeFacesExtNo() const;

  inline const TheMesh& getMesh() const;
  inline const std::vector<std::array<unsigned, 3>>& getBasisComposition() const;

  inline std::vector<FeElement>::const_iterator feElementsCbegin() const;
  inline std::vector<FeElement>::const_iterator feElementsCend() const;
  inline std::vector<FeFaceInt>::const_iterator feFacesIntCbegin() const;
  inline std::vector<FeFaceInt>::const_iterator feFacesIntCend() const;
  inline std::vector<FeFaceExt>::const_iterator feFacesExtCbegin() const;
  inline std::vector<FeFaceExt>::const_iterator feFacesExtCend() const;

  void printElemBasis(std::ostream& out = std::cout) const;
  void printElemBasisDer(std::ostream& out = std::cout) const;
  void printFaceBasis(std::ostream& out = std::cout) const;
  void printFaceBasisDer(std::ostream& out = std::cout) const;

  virtual ~FeSpace() = default;

private:
// Reference to the mesh over which the FeSpace is built
  const TheMesh& Th_;

// Order of polynomials
  unsigned order_;

// Number of degrees of freedom that in 3D is dofNo = (order+1)*(order+2)*(order+3)/(3!)
  unsigned dofNo_;

  std::vector<std::array<unsigned, 3>> basisComposition_;
  std::vector<FeElement> feElements_;
  std::vector<FeFaceInt> feFacesInt_;
  std::vector<FeFaceExt> feFacesExt_;
  const QuadRuleManager::Rule3D& tetraRule_;
  const QuadRuleManager::Rule2D& triaRule_;

// Auxiliary function that computes basisComposition_
  void integerComposition();

// Auxiliary function that fills feElements_, feFacesInt_ and feFacesExt_.
  void initialize();

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline unsigned FeSpace::getOrder() const
{
  return order_;
}

inline unsigned FeSpace::getDofNo() const
{
  return dofNo_;
}

inline const FeElement& FeSpace::getFeElement(unsigned i) const
{
  return feElements_[i];
}

inline const FeFaceInt& FeSpace::getFeFaceInt(unsigned i) const
{
  return feFacesInt_[i];
}

inline const FeFaceExt& FeSpace::getFeFaceExt(unsigned i) const
{
  return feFacesExt_[i];
}

inline unsigned FeSpace::getFeElementsNo() const
{
  return feElements_.size();
}

inline unsigned FeSpace::getFeFacesIntNo() const
{
  return feFacesInt_.size();
}

inline unsigned FeSpace::getFeFacesExtNo() const
{
  return feFacesExt_.size();
}

inline const FeSpace::TheMesh& FeSpace::getMesh() const
{
  return Th_;
}

inline const std::vector<std::array<unsigned, 3>>& FeSpace::getBasisComposition() const
{
  return basisComposition_;
}

inline std::vector<FeElement>::const_iterator FeSpace::feElementsCbegin() const
{
  return feElements_.cbegin();
}

inline std::vector<FeElement>::const_iterator FeSpace::feElementsCend() const
{
  return feElements_.cend();
}

inline std::vector<FeFaceInt>::const_iterator FeSpace::feFacesIntCbegin() const
{
  return feFacesInt_.cbegin();
}

inline std::vector<FeFaceInt>::const_iterator FeSpace::feFacesIntCend() const
{
  return feFacesInt_.cend();
}

inline std::vector<FeFaceExt>::const_iterator FeSpace::feFacesExtCbegin() const
{
  return feFacesExt_.cbegin();
}

inline std::vector<FeFaceExt>::const_iterator FeSpace::feFacesExtCend() const
{
  return feFacesExt_.cend();
}

} // namespace PolyDG

#endif // _FE_SPACE_HPP_
