#ifndef _FE_SPACE_HPP_
#define _FE_SPACE_HPP_

#include "Mesh.hpp"
#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include <vector>

namespace dgfem {

class FeSpace
{
public:
  using theMesh = geom::Mesh;

  FeSpace(theMesh& Th, unsigned order);

  void setOrder(unsigned order);
  unsigned getOrder() const;
  unsigned getDofNo() const;

  virtual ~FeSpace() = default;

private:
  const theMesh& Th_;
  unsigned order_; // ordine dei polinomi
  unsigned dofNo_; // numero di gdl
  std::vector<FeElement> feElements_;
  std::vector<FeFaceInt> feFacesInt_;
  std::vector<FeFaceExt> feFacesExt_;

  void initialize();

};

}

#endif // _FE_SPACE_HPP_
