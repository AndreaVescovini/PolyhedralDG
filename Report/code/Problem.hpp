class Problem
{
public:
  explicit Problem(const FeSpace& Vh);

  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool sym = false);
  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels, bool sym = false);
  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool sym = false);
  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);
  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr, const std::vector<BCLabelType>& bcLabels);

  bool solveLU();
  bool solveCholesky();
  bool solveCG(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
               Real tol = Eigen::NumTraits<Real>::epsilon());
  bool solveBiCGSTAB(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
                     Real tol = Eigen::NumTraits<Real>::epsilon());

  Real computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const;
  Real computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const;

  void exportSolutionVTK(const std::string& fileName, unsigned precision = 8) const;

  // ...

  void finalizeMatrix();

private:
  const FeSpace& Vh_;
  const unsigned dim_;

  Eigen::SparseMatrix<Real> A_;
  Eigen::VectorXd b_;
  Eigen::VectorXd u_;

  std::vector<std::vector<triplet>> triplets_;
  std::vector<bool> sym_;

  // ...
};
