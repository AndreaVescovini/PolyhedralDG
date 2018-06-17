class GradPhiJ : public ExprWrapper<GradPhiJ>
{
public:
  using ReturnType = Eigen::Vector3d;

  const ReturnType& operator()(const FeElement& fe, unsigned /* i */, unsigned j,
                               SizeType t, SizeType p) const
  {
    return fe.getPhiDer(t, p, j);
  }

  const ReturnType& operator()(const FeFaceExt& fe, unsigned /* i */, unsigned j,
                               SizeType p) const
  {
    return fe.getPhiDer(p, j);
  }
};
