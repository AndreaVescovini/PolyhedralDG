class Problem
{
public:
  explicit Problem(const FeSpace& Vh);

  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr,
                    bool sym = false);
  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr,
                         const std::vector<BCLabelType>& bcLab,
                         bool sym = false);
  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr,
                         bool sym = false);
  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);
  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr,
                            const std::vector<BCLabelType>& bcLab);

  void finalizeMatrix();

  // ...

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
