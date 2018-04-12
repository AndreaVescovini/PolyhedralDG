#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "geom.hpp"
#include "ExprWrapper.hpp"
#include <functional>

namespace dgfem
{

class Stiff : public ExprWrapper<Stiff>
{
public:
  Stiff() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~Stiff() = default;
};

class Mass : public ExprWrapper<Mass>
{
public:
  Mass() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~Mass() = default;
};

class GradPhiI : public ExprWrapper<GradPhiI>
{
public:
  GradPhiI() = default;

  const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;

  virtual ~GradPhiI() = default;
};

class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  GradPhiJ() = default;

  const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~GradPhiJ() = default;
};

class PhiI : public ExprWrapper<PhiI>
{
public:
  PhiI() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;

  virtual ~PhiI() = default;
};

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  PhiJ() = default;

  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~PhiJ() = default;
};

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  JumpPhiI() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiI() = default;
};

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  JumpPhiJ() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiJ() = default;
};

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  AverGradPhiI() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiI() = default;
};

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  AverGradPhiJ() = default;

  Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiJ() = default;
};

class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  explicit PenaltyScaling(geom::real sigma = 1.0);

  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~PenaltyScaling() = default;

private:
  geom::real sigma_;
};

class Normal : public ExprWrapper<Normal>
{
public:
  Normal() = default;

  const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, unsigned q) const;
  const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, unsigned j, unsigned q) const;

  virtual ~Normal() = default;
};

class Function : public ExprWrapper<Function>
{
public:
  using fun3real = std::function<geom::real (const Eigen::Vector3d&)>;

  explicit Function(const fun3real& fun);

  geom::real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~Function() = default;

private:
  fun3real fun_;

};

} // namespace dgfem

#endif // _OPERATORS_HPP_
