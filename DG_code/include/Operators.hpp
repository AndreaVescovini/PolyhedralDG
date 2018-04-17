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

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~Stiff() = default;
};

class Mass : public ExprWrapper<Mass>
{
public:
  Mass() = default;

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

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

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  inline geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;

  virtual ~PhiI() = default;
};

class PhiJ : public ExprWrapper<PhiJ>
{
public:
  PhiJ() = default;

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;

  virtual ~PhiJ() = default;
};

class JumpPhiI : public ExprWrapper<JumpPhiI>
{
public:
  JumpPhiI() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiI() = default;
};

class JumpPhiJ : public ExprWrapper<JumpPhiJ>
{
public:
  JumpPhiJ() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline Eigen::Vector3d operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~JumpPhiJ() = default;
};

class AverGradPhiI : public ExprWrapper<AverGradPhiI>
{
public:
  AverGradPhiI() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiI() = default;
};

class AverGradPhiJ : public ExprWrapper<AverGradPhiJ>
{
public:
  AverGradPhiJ() = default;

  inline Eigen::Vector3d operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline const Eigen::Vector3d& operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~AverGradPhiJ() = default;
};

class PenaltyScaling : public ExprWrapper<PenaltyScaling>
{
public:
  explicit PenaltyScaling(geom::real sigma = 1.0);

  inline geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  inline geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~PenaltyScaling() = default;

private:
  geom::real sigma_;
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
  using fun3real = std::function<geom::real (const Eigen::Vector3d&)>;

  explicit Function(const fun3real& fun);

  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const;
  inline geom::real operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const;
  inline geom::real operator()(const FeFaceInt& fe, unsigned i, unsigned j, int side1, int side2, unsigned q) const;
  inline geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned q) const;
  inline geom::real operator()(const FeFaceExt& fe, unsigned i, unsigned j, unsigned q) const;

  virtual ~Function() = default;

private:
  fun3real fun_;

};

//----------------------------------------------------------------------------//
//-------------------------------IMPLEMENTATION-------------------------------//
//----------------------------------------------------------------------------//

inline geom::real Stiff::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhiDer(t, q, i).dot(fe.getPhiDer(t, q, j));
}

inline geom::real Mass::operator()(const FeElement& fe, unsigned i, unsigned j, unsigned t, unsigned q) const
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

inline geom::real PhiI::operator()(const FeElement& fe, unsigned i, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, i);
}

inline geom::real PhiI::operator()(const FeElement& fe, unsigned i, unsigned /* j */, unsigned t, unsigned q) const
{
  return this->operator()(fe, i, t, q);
}

inline geom::real PhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned q) const
{
  return fe.getPhi(q, i);
}

inline geom::real PhiJ::operator()(const FeElement& fe, unsigned /* i */, unsigned j, unsigned t, unsigned q) const
{
  return fe.getPhi(t, q, j);
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, int side1, int /* side2 */, unsigned q) const
{
  return fe.getPhi(side1, q, i) * (1.0 - 2.0 * side1) * fe.getNormal();
}

inline Eigen::Vector3d JumpPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, unsigned q) const
{
  return fe.getPhi(q, i) * fe.getNormal();
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, int /* side1 */, int side2, unsigned q) const
{
  return fe.getPhi(side2, q, j) * (1.0 - 2.0 * side2) * fe.getNormal();
}

inline Eigen::Vector3d JumpPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, unsigned q) const
{
  return fe.getPhi(q, j) * fe.getNormal();
}

inline Eigen::Vector3d AverGradPhiI::operator()(const FeFaceInt& fe, unsigned i, unsigned /* j */, int side1, int /* side2 */, unsigned q) const
{
  return 0.5 * fe.getPhiDer(side1, q, i);
}

inline const Eigen::Vector3d& AverGradPhiI::operator()(const FeFaceExt& fe, unsigned i, unsigned /* j */, unsigned q) const
{
  return fe.getPhiDer(q, i);
}

inline Eigen::Vector3d AverGradPhiJ::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned j, int /* side1 */, int side2, unsigned q) const
{
  return 0.5 * fe.getPhiDer(side2, q, j);
}

inline const Eigen::Vector3d& AverGradPhiJ::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j, unsigned q) const
{
  return fe.getPhiDer(q, j);
}

inline geom::real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* q */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline geom::real PenaltyScaling::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, int /* side1 */, int /* side2 */, unsigned /* q */) const
{
  return sigma_ * fe.getPenaltyParam();
}

inline geom::real PenaltyScaling::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, unsigned /* q */) const
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

inline geom::real Function::operator()(const FeElement& fe, unsigned /* i */, unsigned t, unsigned q) const
{
  return fun_(fe.getQuadPoint(t, q));
}

inline geom::real Function::operator()(const FeElement& fe, unsigned /* i */, unsigned /* j */, unsigned t, unsigned q) const
{
  return fun_(fe.getQuadPoint(t, q));
}

inline geom::real Function::operator()(const FeFaceInt& fe, unsigned /* i */, unsigned /* j */, int /* side1 */, int /* side2 */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

inline geom::real Function::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

inline geom::real Function::operator()(const FeFaceExt& fe, unsigned /* i */, unsigned /* j */, unsigned q) const
{
  return fun_(fe.getQuadPoint(q));
}

} // namespace dgfem

#endif // _OPERATORS_HPP_
