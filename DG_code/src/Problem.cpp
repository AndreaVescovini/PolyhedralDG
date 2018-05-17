#include "Legendre.hpp"
#include "Problem.hpp"
#include "Vertex.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <random>
#include <stdexcept>
#include <unordered_set>

namespace PolyDG
{

Problem::Problem(const FeSpace& Vh)
  : Vh_{Vh}, dim_{static_cast<unsigned>(Vh.getDof() * Vh.getFeElementsNo())},
    A_{dim_, dim_}, b_{Eigen::VectorXd::Zero(dim_)}, u_{Eigen::VectorXd::Zero(dim_)} {}

bool Problem::isSymmetric() const
{
  for(SizeType i = 0; i < sym_.size(); i++)
    if(sym_[i] == false)
      return false;

  return true;
}

void Problem::solveLU()
{
  Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;
  A_.makeCompressed();

  if(this->isSymmetric() == true)
  {
    solver.isSymmetric(true);
    solver.compute(A_.selfadjointView<Eigen::Upper>());
  }
  else
    solver.compute(A_);

  if(solver.info() != Eigen::Success)
    throw std::runtime_error("Numerical issue in the matrix factorization.\n" + solver.lastErrorMessage());

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
    std::cerr << "Numerical issue in the solver.\n" << solver.lastErrorMessage() << std::endl;
}

void Problem::solveChol()
{
  if(this->isSymmetric() == false)
    throw std::domain_error("Error: solveChol() requires a symmetric matrix.");

  Eigen::SimplicialLLT<Eigen::SparseMatrix<Real>, Eigen::Upper> solver;

  A_.makeCompressed();
  solver.compute(A_);
  if(solver.info() != Eigen::Success)
    throw std::runtime_error("Error: Numerical issue in the matrix factorization.");

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
    std::cerr << "Warning: Numerical issue in the solver." << std::endl;
}

void Problem::solveCG(const Eigen::VectorXd& x0, unsigned iterMax, Real tol)
{
  if(this->isSymmetric() == false)
    throw std::domain_error("Error: solveCG() requires a symmetric matrix.");

  Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>, Eigen::Upper> solver;
  solver.setMaxIterations(iterMax);
  solver.setTolerance(tol);

  A_.makeCompressed();
  solver.compute(A_);

  u_ = solver.solveWithGuess(b_, x0);

  if(solver.info() != Eigen::Success)
    std::cerr << "Warning: Conjugate gradient not converged within " << solver.maxIterations() << " iterations." << std::endl;
  else
  {
    // if(verbosity)
    std::cout << "Conjugate gradient converged with " << solver.iterations() << " iterations.\n";
    std::cout << "Estimated error = "<< solver.error() << std::endl;
  }
}

void Problem::solveBiCGSTAB(const Eigen::VectorXd& x0, unsigned iterMax, Real tol)
{
  Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
  solver.setMaxIterations(iterMax);
  solver.setTolerance(tol);

  A_.makeCompressed();

  Eigen::SparseMatrix<Real>  Aselfadj;
  if(this->isSymmetric() == true)
  {
    std::cerr << "Warning: The matrix is symmetric. Consider using the conjugate gradient instead." << std::endl;
    Aselfadj = A_.selfadjointView<Eigen::Upper>();
    solver.compute(Aselfadj);
  }
  else
    solver.compute(A_);

  u_ = solver.solveWithGuess(b_, x0);

  if(solver.info() != Eigen::Success)
    std::cerr << "Warning: BiCGSTAB not converged within " << solver.maxIterations() << " iterations." << std::endl;
  else
  {
    // if(verbosity)
    std::cout << "BiCGSTAB converged with " << solver.iterations() << " iterations.\n";
    std::cout << "Estimated error = "<< solver.error() << std::endl;
  }
}

Real Problem::computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const
{
  Real errSquared = 0.0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    const unsigned indexOffset = it->getElem().getId() * Vh_.getDof();

    for(SizeType t = 0; t < it->getTetrahedraNo(); t++)
      for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
      {
        Real uh = 0.0;

        // Evaluation of the fem function at the quadrature node.
        for(unsigned f = 0; f < Vh_.getDof(); f++)
          uh += u_(f + indexOffset) * it->getPhi(t, p, f);

        const Real difference = uh - uex(it->getQuadPoint(t, p));
        errSquared += difference * difference * it->getWeight(p) * it->getAbsDetJac(t);
      }
  }

  return std::sqrt(errSquared);
}

Real Problem::computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const
{
  Real errSquared = 0.0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    const unsigned indexOffset = it->getElem().getId() * Vh_.getDof();

    for(SizeType t = 0; t < it->getTetrahedraNo(); t++)
      for(SizeType p = 0; p < it->getQuadPointsNo(); p++)
      {
        Eigen::Vector3d uhGrad = Eigen::Vector3d::Zero();

        // Evaluation of the fem function at the quadrature node.
        for(unsigned f = 0; f < Vh_.getDof(); f++)
          uhGrad += u_(f + indexOffset) * it->getPhiDer(t, p, f);

        const Eigen::Vector3d difference = uhGrad - uexGrad(it->getQuadPoint(t, p));
        errSquared += difference.squaredNorm() * it->getWeight(p) * it->getAbsDetJac(t);
      }
  }

  return std::sqrt(errSquared);
}

void Problem::exportSolutionVTK(const std::string& fileName, unsigned precision) const
{
  std::ofstream fout{fileName};

  // Create a vector with random integers in order to distinguish elements.
  std::vector<unsigned> elemValues;
  elemValues.reserve(Vh_.getFeElementsNo());
  for(unsigned i = 0; i < Vh_.getFeElementsNo(); i++)
    elemValues.emplace_back(i);

  std::default_random_engine dre;
  std::shuffle(elemValues.begin(), elemValues.end(), dre);

  // Print the header
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  fout << "  <UnstructuredGrid>\n";

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    const auto& elem = it->getElem();

    // Compute the nodes in the Polyhedron.
    std::unordered_set<std::reference_wrapper<const Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> nodesSet;
    nodesSet.reserve(elem.getTetrahedraNo() + 3);

    for(SizeType i = 0; i < elem.getTetrahedraNo(); i++)
      for(SizeType j = 0; j < 4; j++)
        nodesSet.emplace(elem.getTetra(i).getVertex(j));

    const std::vector<std::reference_wrapper<const Vertex>> nodes(nodesSet.cbegin(), nodesSet.cend());

    // Compute the solution at the nodes.
    std::vector<Real> uNodes;
    uNodes.reserve(nodes.size());
    for(auto itNod = nodes.cbegin(); itNod != nodes.cend(); itNod++)
      uNodes.emplace_back(evalSolution(itNod->get().getX(), itNod->get().getY(), itNod->get().getZ(), *it));

    fout << "    <Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << elem.getTetrahedraNo() << "\">\n";

    fout << std::setprecision(precision) << std::scientific;

    // Print the nodes coordinates.
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n         ";
    for(auto itNod = nodes.cbegin(); itNod != nodes.cend(); itNod++)
      fout << ' ' << itNod->get().getX() << ' ' << itNod->get().getY() << ' ' << itNod->get().getZ();
    fout << "\n        </DataArray>\n";
    fout << "      </Points>\n";

    // Print the cells (tetrahedra) connectivity and type.
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n         ";
    for(SizeType i = 0; i < elem.getTetrahedraNo(); i++)
      for(SizeType j = 0; j < 4; j++)
        fout << ' ' << (std::find(nodes.cbegin(), nodes.cend(), elem.getTetra(i).getVertex(j)) - nodes.cbegin());
    fout << "\n        </DataArray>\n";

    fout << "         <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n         ";
    for(unsigned offset = 4; offset <= elem.getTetrahedraNo() * 4; offset += 4)
      fout << ' ' << offset;
    fout << "\n        </DataArray>\n";

    fout << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < elem.getTetrahedraNo(); i++)
      fout << " 10";
    fout << "\n        </DataArray>\n";
    fout << "      </Cells>\n";

    // Print the values of the solution.
    fout << "      <PointData Scalars=\"Solution\">\n";
    fout << "        <DataArray type=\"Float64\" Name=\"Solution\" format=\"ascii\">\n         ";
    for(SizeType i = 0; i < uNodes.size(); i++)
      fout << ' ' << uNodes[i];
    fout << "\n        </DataArray>\n";
    fout << "      </PointData>\n";

    // Print a value for the Polyhedron.
    fout << "      <CellData Scalars=\"Mesh\">\n";
    fout << "        <DataArray type=\"UInt32\" Name=\"Mesh\" format=\"ascii\">\n         ";
    for(SizeType i = 0; i < elem.getTetrahedraNo(); i++)
      fout << ' ' << elemValues[elem.getId()];
    fout << "\n        </DataArray>\n";
    fout << "      </CellData>\n";

    fout << "    </Piece>\n";
  }

  fout << "  </UnstructuredGrid>\n";
  fout << "</VTKFile>" << std::endl;

  fout.close();
}

void Problem::printInfo(std::ostream& out) const
{
  out << "-------------------- PROBLEM INFO --------------------" << '\n';
  out << "Symmetric:  " << (this->isSymmetric() == true ? "Yes" : "No") << '\n';
  out << "Total degrees of freedom: " << dim_ <<'\n';
  if(A_.nonZeros() != 0)
    out << "Actual non-zeros in the matrix: " << A_.nonZeros() << '\n';
  out << "------------------------------------------------------" << std::endl;
}

Real Problem::evalSolution(Real x, Real y, Real z, const FeElement& el) const
{
  Real result = 0.0;
  const unsigned indexOffset = el.getElem().getId() * Vh_.getDof();
  const Eigen::Vector3d hb = el.getElem().getBoundingBox().sizes() / 2;
  const Eigen::Vector3d mb = el.getElem().getBoundingBox().center();

  const auto& basisComposition = Vh_.getBasisComposition();

  for(unsigned i = 0; i < Vh_.getDof(); i++)
  {
    const Real valx = legendre(basisComposition[i][0], (x - mb(0)) / hb(0)) / std::sqrt(hb(0));
    const Real valy = legendre(basisComposition[i][1], (y - mb(1)) / hb(1)) / std::sqrt(hb(1));
    const Real valz = legendre(basisComposition[i][2], (z - mb(2)) / hb(2)) / std::sqrt(hb(2));

    result += u_(indexOffset + i) * valx  * valy * valz;
  }

  return result;
}

void Problem::finalizeMatrix()
{
  std::vector<triplet> concatVec;

  SizeType dimTot = 0;
  for(SizeType i = 0; i < triplets_.size(); i++)
    dimTot += (sym_[i] == true ? 2 : 1) * triplets_[i].size();

  concatVec.reserve(dimTot);

  for(SizeType i = 0; i < triplets_.size(); i++)
  {
    if(this->isSymmetric() == false && sym_[i] == true)
      for(const auto& t : triplets_[i])
        if(t.row() != t.col())
          concatVec.emplace_back(t.col(), t.row(), t.value());

    concatVec.insert(concatVec.cend(),
                     std::make_move_iterator(triplets_[i].cbegin()),
                     std::make_move_iterator(triplets_[i].cend()));
    triplets_[i].clear();
  }

  A_.setFromTriplets(concatVec.cbegin(), concatVec.cend());
  A_.prune(A_.coeff(0,0));

  triplets_.clear();
}

void Problem::clearMatrix()
{
  A_.setZero();
  sym_.clear();
}

void Problem::clearRhs()
{
  b_ = Eigen::VectorXd::Zero(dim_);
}

} // namespace PolyDG
