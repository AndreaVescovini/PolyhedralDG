#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "geom.hpp"
#include "ExprWrapper.hpp"

namespace dgfem
{

class Stiff : public ExprWrapper<Stiff>
{
public:
  Stiff() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;
  virtual ~Stiff() = default;

};

class Mass : public ExprWrapper<Mass>
{
public:
  Mass() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j,
                        unsigned t, unsigned q) const;
  virtual ~Mass() = default;

};

class GradPhi : public ExprWrapper<GradPhi>
{
public:
  GradPhi() = default;

  Eigen::Vector3d operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;

  virtual ~GradPhi() = default;
};

class Phi : public ExprWrapper<Phi>
{
public:
  Phi() = default;

  geom::real operator()(const FeElement&, unsigned i, unsigned t, unsigned q) const;

  virtual ~Phi() = default;
};

class PhiI : public ExprWrapper<PhiI>
{
public:
  PhiI() = default;

  geom::real operator()(const FeElement&, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~PhiI() = default;
};

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  PhiJ() = default;

  geom::real operator()(const FeElement&, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~PhiJ() = default;
};

class JumpInt : public ExprWrapper<JumpInt>
{
public:
  JumpInt() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned side, unsigned q) const;

  virtual ~JumpInt() = default;
};

}

#endif // _OPERATORS_HPP_
