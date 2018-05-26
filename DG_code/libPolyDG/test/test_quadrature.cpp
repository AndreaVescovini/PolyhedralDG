// Test for the exactness of the quadrature rules integrating polynomials

#include "PolyDG.hpp"
#include "QuadRuleManager.hpp"
#include "Utilities.hpp"
#include "Watch.hpp"

#include <Eigen/Core>

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

int main()
{
  using PolyDG::QuadRuleManager;
  using PolyDG::Real;
  using PolyDG::SizeType;
  using Utilities::pow;

  Utilities::Watch ch;
  ch.start();

  const Real tol = 1e-14;

  const auto& qm = QuadRuleManager::instance();

  std::cout << "Testing weights of rules over tetrahedra..." << '\n';
  for(auto it = qm.tetraCbegin(); it != qm.tetraCend(); it++)
    std::cout << "Rule " << it->getDoe() << (it->checkRuleWeights(tol) ? " is ok." : " is wrong.") << '\n';

  std::cout << "\nTesting weights of rules over triangles..." << '\n';
  for(auto it = qm.triaCbegin(); it != qm.triaCend(); it++)
    std::cout << "Rule " << it->getDoe() << (it->checkRuleWeights(tol) ? " is ok." : " is wrong.") << '\n';

  std::vector<std::function<PolyDG::Real (const Eigen::Vector3d&)>> funTetra;
  funTetra.reserve(qm.getTetraRuleNo());
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 24.0 * x(0); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 120.0 * x(0) * x(1); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 720.0 * x(0) * x(1) * x(2); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 2520.0 * x(0) * x(0) * x(1) * x(2); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 6720.0 * pow(x(0), 3) * x(1) * x(2); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 30240.0 * pow(x(0), 3) * x(1) * x(1) * x(2); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 151200.0 * pow(x(0) * x(1) * x(2), 2) * x(0); });
  funTetra.emplace_back([](const Eigen::Vector3d& x) { return 554400.0 * pow(x(0) * x(1), 3) * x(2) * x(2); });

  std::cout << std::setprecision(-std::log10(tol)) << std::fixed;

  std::cout << "\nTesting exactness of rules over tetrahedra..." << '\n';
  for(auto it = qm.tetraCbegin(); it != qm.tetraCend(); it++)
  {
    Real sum = 0.0;
    for(SizeType i = 0; i < it->getPointsNo(); i++)
      sum += funTetra[it->getDoe() - 1](it->getPoint(i)) * it->getWeight(i);

    std::cout << "Rule " << it->getDoe() << ", expected value = 1, computed value = " << sum  <<  (std::abs(sum - 1.0) < tol ? ",  ok." : ",  wrong.") << '\n';
  }

  std::vector<std::function<PolyDG::Real (const Eigen::Vector2d&)>> funTria;
  funTria.reserve(qm.getTetraRuleNo());
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 6.0 * x(0); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 24.0 * x(0) * x(1); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 60.0 * x(0) * x(0) * x(1); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 180.0 * pow(x(0) * x(1), 2); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 420.0 * pow(x(0) * x(1), 2) * x(0); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 840.0 * pow(x(0), 4) * x(1) * x(1); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 2520.0 * pow(x(0) * x(1), 3) * x(0); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 6300.0 * pow(x(0) * x(1), 4); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 13860.0 * pow(x(0) * x(1), 4) * x(0); });
  funTria.emplace_back([](const Eigen::Vector2d& x) { return 33264.0 * pow(x(0) * x(1), 5); });

  std::cout << "\nTesting exactness of rules over triangles..." << '\n';
  for(auto it = qm.triaCbegin(); it != qm.triaCend(); it++)
  {
    Real sum = 0.0;
    for(SizeType i = 0; i < it->getPointsNo(); i++)
      sum += funTria[it->getDoe() - 1](it->getPoint(i)) * it->getWeight(i);

      std::cout << "Rule " << it->getDoe() << ", expected value = 1, computed value = " << sum << (std::abs(sum - 1.0) < tol ? ",  ok." : ",  wrong.") << '\n';
  }

  std::cout << "\nTest finished." << std::endl;

  ch.stop();
  std::cout << ch << std::endl;

  return 0;
}
