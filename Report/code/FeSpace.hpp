class FeSpace
{
public:
  FeSpace(Mesh& Th, unsigned degree, unsigned doeQuad3D, unsigned doeQuad2D);
  FeSpace(Mesh& Th, unsigned degree);

  // Getter methods

  // ...

  // Alias for a random access const iterator over FeElement, FeFaceExt or
  // FeFaceInt
  template <typename T>
  using ConstIter = typename std::vector<T>::const_iterator;

  inline ConstIter<FeElement> feElementsCbegin() const;
  inline ConstIter<FeElement> feElementsCend() const;
  // The same for FeFaceExt and FeFaceInt...

private:
  const Mesh& Th_;
  const unsigned degree_;
  const unsigned dof_;

  // Possible degrees of the monomials that multiplied togheter give polynomials
  // of degree less or equal to degree_
  std::vector<std::array<unsigned, 3>> basisComposition_;

  std::vector<FeElement> feElements_;
  std::vector<FeFaceExt> feFacesExt_;
  std::vector<FeFaceInt> feFacesInt_;

  const QuadRule3D& tetraRule_;
  const QuadRule2D& triaRule_;

  void integerComposition();
  void initialize();
};
