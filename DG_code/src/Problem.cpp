#include "Problem.hpp"
#include "Vertex.hpp"
#include "LegendreBasis.hpp"

#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

#include <cmath>
#include <fstream>
#include <unordered_set>
#include <algorithm>
#include <iterator>

namespace PolyDG
{

Problem::Problem(const FeSpace& Vh)
  : Vh_{Vh}, dim_{Vh.getFeElementsNo() * Vh.getDofNo()}, A_{dim_, dim_},
    b_{Eigen::VectorXd::Zero(dim_)}, u_{Eigen::VectorXd::Zero(dim_)},
    sym_{true, true, true} {}

void Problem::solveLU()
{
  Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;

  if(this->isSymmetric() == true)
  {
    solver.isSymmetric(true);
    solver.compute(A_.selfadjointView<Eigen::Upper>());
  }
  else
    solver.compute(A_);

  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue\n" << solver.lastErrorMessage() << std::endl;
    return;
  }

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue\n" << solver.lastErrorMessage() << std::endl;
    return;
  }
}

void Problem::solveChol()
{
  if(this->isSymmetric() == false)
  {
    std::cerr << "solveChol() requires a symmetric matrix" << std::endl;
    return;
  }

  // Eigen::SimplicialLLT<Eigen::SparseMatrix<Real>, Eigen::Upper> solver;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>, Eigen::Upper> solver;

  A_.makeCompressed();
  solver.compute(A_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue " << solver.info() << std::endl;
    return;
  }

  u_ = solver.solve(b_);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Numerical Issue" << std::endl;
    return;
  }
}

void Problem::solveCG(const Eigen::VectorXd& x0, unsigned iterMax, Real tol)
{
  if(this->isSymmetric() == false)
  {
    std::cerr << "solveCG() requires a symmetric matrix" << std::endl;
    return;
  }

  Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>, Eigen::Upper> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>, Eigen::Lower|Upper> solver;

  solver.setMaxIterations(iterMax);
  solver.setTolerance(tol);

  solver.compute(A_);

  u_ = solver.solveWithGuess(b_, x0);
  if(solver.info() != Eigen::Success)
  {
    std::cerr << "Not converged" << std::endl;
    return;
  }

  // if(verbosity)
  std::cout << "Converged with " << solver.iterations() << " iterations.\n";
  std::cout << "Estimated error "<< solver.error() << std::endl;
}

Real Problem::computeErrorL2(const std::function<Real (const Eigen::Vector3d&)>& uex) const
{
  Real errSquared = 0.0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = it->getElem().getId() * Vh_.getDofNo();

    for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
      for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
      {
        Real uh = 0.0;

        // Evaluation of the fem function at the quadrature node
        for(unsigned i = 0; i < Vh_.getDofNo(); i++)
          uh += u_(i + indexOffset) * it->getPhi(t, q, i);

        Real difference = uh - uex(it->getQuadPoint(t, q));
        errSquared += difference * difference * it->getWeight(q) * it->getAbsDetJac(t);
      }
  }

  return std::sqrt(errSquared);
}

Real Problem::computeErrorH10(const std::function<Eigen::Vector3d (const Eigen::Vector3d&)>& uexGrad) const
{
  Real errSquared = 0.0;

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    unsigned indexOffset = it->getElem().getId() * Vh_.getDofNo();

    for(unsigned t = 0; t < it->getTetrahedraNo(); t++)
      for(unsigned q = 0; q < it->getQuadPointsNo(); q++)
      {
        Eigen::Vector3d uhGrad = Eigen::Vector3d::Zero();

        // Evaluation of the fem function at the quadrature node
        for(unsigned i = 0; i < Vh_.getDofNo(); i++)
          uhGrad += u_(i + indexOffset) * it->getPhiDer(t, q, i);

        Eigen::Vector3d difference = uhGrad - uexGrad(it->getQuadPoint(t, q));
        errSquared += difference.squaredNorm() * it->getWeight(q) * it->getAbsDetJac(t);
      }
  }

  return std::sqrt(errSquared);
}

void Problem::exportSolutionVTK(const std::string& fileName) const
{
  std::ofstream fout{fileName};

  // Print the header
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  fout << "  <UnstructuredGrid>\n";

  for(auto it = Vh_.feElementsCbegin(); it != Vh_.feElementsCend(); it++)
  {
    auto elem = it->getElem();

  // Compute the nodes in the Polyhedron
    std::unordered_set<std::reference_wrapper<const Vertex>, std::hash<Vertex>, std::equal_to<Vertex>> nodesSet;
    nodesSet.reserve(elem.getTetrahedraNo() + 3);

    for(unsigned i = 0; i < elem.getTetrahedraNo(); i++)
      for(unsigned j = 0; j < 4; j++)
        nodesSet.emplace(elem.getTetra(i).getVertex(j));

    std::vector<std::reference_wrapper<const Vertex>> nodes(nodesSet.cbegin(), nodesSet.cend());

  // Compute the solution at the nodes
    std::vector<Real> uNodes(nodes.size());
    uNodes.reserve(nodes.size());
    for(auto itNod = nodes.cbegin(); itNod != nodes.cend(); itNod++)
      uNodes.emplace_back(evalSolution(itNod->get().getX(), itNod->get().getY(), itNod->get().getZ(), *it));

    fout << "    <Piece NumberOfPoints=\"" << nodes.size() << "\" NumberOfCells=\"" << elem.getTetrahedraNo() << "\">\n";

  // Print the nodes coordinates
    fout << "      <Points>\n";
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n         ";
    for(auto itNod = nodes.cbegin(); itNod != nodes.cend(); itNod++)
      fout << ' ' << itNod->get().getX() << ' ' << itNod->get().getY() << ' ' << itNod->get().getZ();
    fout << "\n        </DataArray>\n";
    fout << "      </Points>\n";

  // Print the cells (tetrahedra) connectivity and type
    fout << "      <Cells>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < elem.getTetrahedraNo(); i++)
      for(unsigned j = 0; j < 4; j++)
      {
        unsigned vertexIndex = std::find(nodes.cbegin(), nodes.cend(), elem.getTetra(i).getVertex(j)) - nodes.cbegin();
        fout << ' ' << vertexIndex;
      }
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

  // Print the values of the solution
    fout << "      <PointData Scalars=\"Solution\">\n";
    fout << "        <DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < uNodes.size(); i++)
      fout << ' ' << uNodes[i];
    fout << "\n        </DataArray>\n";
    fout << "      </PointData>\n";

  // Print a value for the Polyhedron
    fout << "      <CellData Scalars=\"Mesh\">\n";
    fout << "        <DataArray type=\"Int32\" Name=\"Mesh\" format=\"ascii\">\n         ";
    for(unsigned i = 0; i < elem.getTetrahedraNo(); i++)
      fout << ' ' << elem.getId();
    fout << "\n        </DataArray>\n";
    fout << "      </CellData>\n";

    fout << "    </Piece>\n";
  }

  fout << "  </UnstructuredGrid>\n";
  fout << "</VTKFile>" << std::endl;

  fout.close();
}

Real Problem::evalSolution(Real x, Real y, Real z, const FeElement& el) const
{
  Real result = 0.0;
  // unsigned indexOffset = el.getElem().getId() * Vh_.getDofNo();
  // Eigen::Vector3d hb = el.getElem().getBoundingBox().sizes() / 2;
  // Eigen::Vector3d mb = el.getElem().getBoundingBox().center();
  //
  // auto basisComposition = Vh_.getBasisComposition();

  // for(unsigned i = 0; i < Vh_.getDofNo(); i++)
  // {
  //   Real valx = LegendreBasis<basisComposition[i][0]>((x - mb(0)) / hb(0)) / std::sqrt(hb(0));
  //   Real valy = LegendreBasis<basisComposition[i][1]>((y - mb(1)) / hb(1)) / std::sqrt(hb(1));
  //   Real valz = LegendreBasis<basisComposition[i][2]>((z - mb(2)) / hb(2)) / std::sqrt(hb(2));
  //
  //   result += u_(indexOffset + i) * valx  * valy * valz;
  // }

  return result;
}

void Problem::finalizeMatrix()
{
  std::vector<triplet> concatVec;

  concatVec.reserve((sym_[0] == true ? 2 : 1) * triplets_[0].size() +
                    (sym_[1] == true ? 2 : 1) * triplets_[1].size() +
                    (sym_[2] == true ? 2 : 1) * triplets_[2].size());

  for(unsigned i = 0; i < 3; i++)
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
}

} // namespace PolyDG
