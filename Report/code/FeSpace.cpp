FeSpace::FeSpace(Mesh& Th, unsigned degree, unsigned doeQuad3D, unsigned doeQuad2D)
  : Th_{Th}, degree_{degree}, dof_{(degree + 1) * (degree + 2) * (degree + 3) / 6},
    tetraRule_{QuadRuleManager::instance().getTetraRule(doeQuad3D)},
    triaRule_ {QuadRuleManager::instance().getTriaRule(doeQuad2D)}
  {
    integerComposition();
    initialize();
  }

FeSpace::FeSpace(Mesh& Th, unsigned degree)
  : FeSpace(Th, degree, 2 * (degree - 1), 2 * degree) {}

// ...

void FeSpace::initialize()
{
  feElements_.reserve(Th_.getPolyhedraNo());
  for(SizeType i = 0; i < Th_.getPolyhedraNo(); i++)
    feElements_.emplace_back(Th_.getPolyhedron(i), dof_, basisComposition_, tetraRule_);

  // The same for feFacesExt_ and feFacesInt_
}
