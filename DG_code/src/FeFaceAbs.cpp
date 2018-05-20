/*!
    @file   FeFaceAbs.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class FeFaceAbs
*/

#include "FeFaceAbs.hpp"

namespace PolyDG
{

FeFaceAbs::FeFaceAbs(const FaceAbs& face, unsigned dof,
                     const std::vector<std::array<unsigned, 3>>& basisComposition,
                     const QuadRule2D& triaRule)
  : face_{face}, dof_{dof}, basisComposition_{basisComposition},
    triaRule_{triaRule} {}

} // namespace PolyDG
