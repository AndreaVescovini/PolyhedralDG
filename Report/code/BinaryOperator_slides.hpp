template <typename LO, typename RO, typename OP>
class BinaryOperator : public ExprWrapper<BinaryOperator<LO, RO, OP>> {
public:
  // Alias of return type of the binary operator
  using ReturnType = decltype(OP()(std::declval<typename LO::ReturnType>(),
                                   std::declval<typename RO::ReturnType>()));

  ReturnType operator()(const FeElement& fe, unsigned i, unsigned j, SizeType t,
                        SizeType p) const {
    return OP()(lo_(fe, i, j, t, p), ro_(fe, i, j, t, p));
  }
  // ... [Overload for FeFaceExt and FeFaceInt]

private:
  const LO& lo_;
  const RO& ro_;
};
