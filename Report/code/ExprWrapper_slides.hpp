template <typename E> class ExprWrapper {
public:
  // ...
  operator E&() { return *static_cast<E*>(this); }
  E& asDerived() { return *static_cast<E*>(this); }
};
