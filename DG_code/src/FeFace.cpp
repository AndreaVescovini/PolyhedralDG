#include "FeFace.hpp"

namespace PolyDG
{

FeFace::FeFace(unsigned dof,
               const std::vector<std::array<unsigned, 3>>& basisComposition,
               const QuadRule<Eigen::Vector2d>& triaRule)
  : dof_{dof}, basisComposition_{basisComposition}, triaRule_{triaRule} {}

} // namespace PolyDG
