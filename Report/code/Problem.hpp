class Problem
{
public:
  explicit Problem(const FeSpace& Vh);

  template <typename T>
  void integrateVol(const ExprWrapper<T>& expr, bool sym = false);
  template <typename T>
  void integrateFacesExt(const ExprWrapper<T>& expr,
                         const std::vector<BCLabelType>& bcLabels, bool sym=false);
  template <typename T>
  void integrateFacesInt(const ExprWrapper<T>& expr, bool sym = false);
  template <typename T>
  void integrateVolRhs(const ExprWrapper<T>& expr);
  template <typename T>
  void integrateFacesExtRhs(const ExprWrapper<T>& expr,
                            const std::vector<BCLabelType>& bcLabels);

  bool solveLU();
  bool solveCholesky();
  bool solveCG(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
               Real tol = Eigen::NumTraits<Real>::epsilon());
  bool solveBiCGSTAB(const Eigen::VectorXd& x0, unsigned iterMax = 10000,
                     Real tol = Eigen::NumTraits<Real>::epsilon());

  Real computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const;
  Real computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const;

  void exportSolutionVTK(const std::string& fileName, unsigned precision=8) const;

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

template <typename T>
void Problem::integrateVol(const ExprWrapper<T>& expr, bool sym)
{
 // ...

 for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
 {
   const unsigned indexOffset = it->getElem().getId() * Vh_.getDof();
   for(unsigned j = 0; j < Vh_.getDof(); j++)
     for(unsigned i = 0; i < (sym == true ? j + 1 : Vh_.getDof()); i++)
     {
       Real sum = 0.0;
       for(SizeType t = 0; t < it->getTetrahedraNo(); t++)
         for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
           sum += exprDerived(*it,i,j,t,p) * it->getWeight(p)*it->getAbsDetJac(t);

         triplets_.back().emplace_back(i + indexOffset, j + indexOffset, sum);
     }
 }
}
