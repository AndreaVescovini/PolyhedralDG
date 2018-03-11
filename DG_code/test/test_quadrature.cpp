#include <iostream>
#include "QuadRuleManager.hpp"

using namespace dgfem;

int main()
{
  std::cout << "Testing rules over tetrahedra..." << '\n';
  for(auto it = QuadRuleManager::tetraCbegin(); it != QuadRuleManager::tetraCend(); it++)
  {
    if(it->checkRule() == true)
      std::cout << "Rule " << it->getDoe() << " is ok." << '\n';
    else
      std::cout << "Rule " << it->getDoe() << " is wrong." << '\n';
  }

  std::cout << "\nTesting rules over triangles..." << '\n';
  for(auto it = QuadRuleManager::triaCbegin(); it != QuadRuleManager::triaCend(); it++)
  {
    if(it->checkRule() == true)
      std::cout << "Rule " << it->getDoe() << " is ok." << '\n';
    else
      std::cout << "Rule " << it->getDoe() << " is wrong." << '\n';
  }
  std::cout << "Test finished" << std::endl;


  return 0;
}
