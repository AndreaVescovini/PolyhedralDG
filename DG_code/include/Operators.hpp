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

class Stiff : public ExprWrapper<Stiff>
{
public:
  Stiff() = default;
  Stiff(const Stiff&) = default;
  Stiff(Stiff&&) = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  virtual ~Stiff() = default;
};

class Mass : public ExprWrapper<Mass>
{
public:
  Mass() = default;
  Mass(const Mass&) = default;
  Mass(Mass&&) = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  virtual ~Mass() = default;
};

class GradPhiI : public ExprWrapper<GradPhiI>
{
public:
  GradPhiI() = default;
  GradPhiI(const GradPhiI&) = default;
  GradPhiI(GradPhiI&&) = default;

  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  virtual ~GradPhiI() = default;
};

class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  GradPhiJ() = default;
  GradPhiJ(const GradPhiJ&) = default;
  GradPhiJ(GradPhiJ&&) = default;

  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  virtual ~GradPhiJ() = default;
};

class PhiI : public ExprWrapper<PhiI>
{
public:
  PhiI() = default;
  PhiI(const PhiI&) = default;
  PhiI(PhiI&&) = default;

  inline Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;

  virtual ~PhiI() = default;
};

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  PhiJ() = default;
  PhiJ(const PhiJ&) = default;
  PhiJ(PhiJ&&) = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;

  virtual ~PhiJ() = default;
};

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  JumpPhiI() = default;
  JumpPhiI(const JumpPhiI&) = default;
  JumpPhiI(JumpPhiI&&) = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~JumpPhiI() = default;
};

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  JumpPhiJ() = default;
  JumpPhiJ(const JumpPhiJ&) = default;
  JumpPhiJ(JumpPhiJ&&) = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~JumpPhiJ() = default;
};

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  AverGradPhiI() = default;
  AverGradPhiI(const AverGradPhiI&) = default;
  AverGradPhiI(AverGradPhiI&&) = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~AverGradPhiI() = default;
};

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  AverGradPhiJ() = default;
  AverGradPhiJ(const AverGradPhiJ&) = default;
  AverGradPhiJ(AverGradPhiJ&&) = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~AverGradPhiJ() = default;
};

class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  explicit PenaltyScaling(Real sigma = 1.0);
  PenaltyScaling(const PenaltyScaling&) = default;
  PenaltyScaling(PenaltyScaling&&) = default;

  inline Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  inline Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~PenaltyScaling() = default;

private:
  Real sigma_;
};

class Normal : public ExprWrapper<Normal>
{
public:
  Normal() = default;
  Normal(const Normal&) = default;
  Normal(Normal&&) = default;

  inline const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, SizeType p) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, unsigned j, SizeType p) const;

  virtual ~Normal() = default;
};

class Function : public ExprWrapper<Function>
{
public:
  using fun3real = std::function<Real (const Eigen::Vector3d&)>;

  explicit Function(const fun3real& fun);
  Function(const Function&) = default;
  Function(Function&&) = default;

  inline Real operator()(const FeElement& fe, unsigned i, SizeType t, SizeType p) const;
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t, SizeType p) const;
  inline Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, SizeType p) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, SizeType p) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, SizeType p) const;

  virtual ~Function() = default;

private:
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
  return this->operator()(fe, i, t, p);
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
