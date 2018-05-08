#include "FeFace.hpp"

namespace PolyDG
{

FeFace::FeFace(const FaceAbs& face, unsigned dof,
               const std::vector<std::array<unsigned, 3>>& basisComposition,
               const QuadRule<Eigen::Vector2d>& triaRule)
  : face_{face}, dof_{dof}, basisComposition_{basisComposition},
    triaRule_{triaRule} {}

} // namespace PolyDG
