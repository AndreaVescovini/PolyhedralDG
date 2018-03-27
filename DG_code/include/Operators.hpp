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

class JumpIntPhi : public ExprWrapper<JumpIntPhi>
{
public:
  JumpIntPhi() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, int side, unsigned q) const;

  virtual ~JumpIntPhi() = default;
};

class AverIntGradPhi : public ExprWrapper<AverIntGradPhi>
{
public:
  AverIntGradPhi() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, int side, unsigned q) const;

  virtual ~AverIntGradPhi() = default;
};

class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  explicit PenaltyScaling(geom::real sigma = 1.0);

  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;

  virtual ~PenaltyScaling() = default;

private:
  geom::real sigma_;
};

}

#endif // _OPERATORS_HPP_
