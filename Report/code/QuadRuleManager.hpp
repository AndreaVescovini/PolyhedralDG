class QuadRuleManager
{
public:
  static QuadRuleManager& instance();

  // ...

  void setTetraRule(const QuadRule3D& rule);
  void setTetraRule(QuadRule3D&& rule);
  void setTriaRule(const QuadRule2D& rule);
  void setTriaRule(QuadRule2D&& rule);

private:
  QuadRuleManager();

  std::set<QuadRule3D, std::less<QuadRule3D>> tetraRules_;
  std::set<QuadRule2D, std::less<QuadRule2D>> triaRules_;

  std::array<Eigen::Matrix3d, 4> faceMaps_;
};
