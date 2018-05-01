#ifndef _OPERATORS_HPP_
#define _OPERATORS_HPP_

#include "PolyDG.hpp"
#include "FeElement.hpp"
#include "FeFaceInt.hpp"
#include "FeFaceExt.hpp"
#include "ExprWrapper.hpp"

#include <functional>

namespace PolyDG
{

class Stiff : public ExprWrapper<Stiff>
{
public:
  Stiff() = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~Stiff() = default;
};

class Mass : public ExprWrapper<Mass>
{
public:
  Mass() = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~Mass() = default;
};

class GradPhiI : public ExprWrapper<GradPhiI>
{
public:
  GradPhiI() = default;

  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;

  virtual ~GradPhiI() = default;
};

class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  GradPhiJ() = default;

  inline const Eigen::Vector3d& operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~GradPhiJ() = default;
};

class PhiI : public ExprWrapper<PhiI>
{
public:
  PhiI() = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;

  virtual ~PhiI() = default;
};

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  PhiJ() = default;

  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~PhiJ() = default;
};

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  JumpPhiI() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiI() = default;
};

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  JumpPhiJ() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiJ() = default;
};

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  AverGradPhiI() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiI() = default;
};

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  AverGradPhiJ() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiJ() = default;
};

class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  explicit PenaltyScaling(Real sigma = 1.0);

  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  inline Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~PenaltyScaling() = default;

private:
  Real sigma_;
};

class Normal : public ExprWrapper<Normal>
{
public:
  Normal() = default;

  inline const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt&, unsigned i, unsigned j, unsigned q) const;

  virtual ~Normal() = default;
};

class Function : public ExprWrapper<Function>
{
public:
  using fun3real = std::function<Real (const Eigen::Vector3d&)>;

  explicit Function(const fun3real& fun);

  inline Real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  inline Real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  inline Real operator()(const FeFaceInt& fe, unsigned i, unsigned j, SideType si, SideType sj, unsigned q) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  inline Real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~Function() = default;

private:
  fun3real fun_;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline Real Stiff::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i).dot(fe.getPhiDer(t, q, j));
}

inline Real Mass::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i) * fe.getPhi(t, q, j);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeElement& fe, unsigned i, unsigned /* j */, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i);
}

inline const Eigen::Vector3d& GradPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return fe.getPhiDer(q, i);
}

inline const Eigen::Vector3d& GradPhiJ::operator()(const FeElement& fe, unsigned /* i */, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, j);
}

inline Real PhiI::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i);
}

inline Real PhiI::operator()(const FeElement& fe, unsigned i, unsigned /* j */, unsigned t, unsigned q) const
{
  return this->operator()(fe, i, t, q);
}

inline Real PhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return fe.getPhi(q, i);
}

inline Real PhiJ::operator()(const FeElement& fe, unsigned /* i */, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, j);
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, unsigned q) const
{
  return fe.getPhi(si, q, i) * fe.getNormal() * (si == Out ? 1 : -1);
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, unsigned q) const
{
  return fe.getPhi(q, i) * fe.getNormal();
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, unsigned q) const
{
  return fe.getPhi(sj, q, j) * fe.getNormal() * (sj == Out ? 1 : -1);
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, unsigned q) const
{
  return fe.getPhi(q, j) * fe.getNormal();
}

inline Eigen::Vector3d AverGradPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, SideType si, SideType /* sj */, unsigned q) const
{
  return 0.5 * fe.getPhiDer(si, q, i);
}

inline const Eigen::Vector3d& AverGradPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, unsigned q) const
{
  return fe.getPhiDer(q, i);
}

inline Eigen::Vector3d AverGradPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, SideType /* si */, SideType sj, unsigned q) const
{
  return 0.5 * fe.getPhiDer(sj, q, j);
}

inline const Eigen::Vector3d& AverGradPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, unsigned q) const
{
  return fe.getPhiDer(q, j);
}

inline Real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* q */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline Real PenaltyScaling::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, unsigned /* q */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline Real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, unsigned /* q */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline const Eigen::Vector3d& Normal::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* q */) const
{
  return fe.getNormal();
}

inline const Eigen::Vector3d& Normal::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, unsigned /* q */) const
{
  return fe.getNormal();
}

inline Real Function::operator()(const FeElement& fe, unsigned /* i */, unsigned t, unsigned q) const
{
  return fun_(fe.getQuadPoint(t, q));
}

inline Real Function::operator()(const FeElement& fe, unsigned /* i */, unsigned /* j */, unsigned t, unsigned q) const
{
  return fun_(fe.getQuadPoint(t, q));
}

inline Real Function::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, SideType /* si */, SideType /* sj */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

inline Real Function::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

inline Real Function::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

} // namespace PolyDG

#endif // _OPERATORS_HPP_
