#include "QuadRuleManager.hpp"

#include <algorithm>

namespace PolyDG
{

QuadRuleManager& QuadRuleManager::instance()
{
  static QuadRuleManager qrm;
  return qrm;
}

const QuadRule3D& QuadRuleManager::getTetraRule(unsigned doe) const
{
  ConstIter<QuadRule3D> iter = std::lower_bound(tetraRules_.cbegin(), tetraRules_.cend(), doe, compDoe<QuadRule3D>);
  return iter != tetraRules_.cend() ? *iter : *(--iter);
}

const QuadRule2D& QuadRuleManager::getTriaRule(unsigned doe) const
{
  ConstIter<QuadRule2D> iter = std::lower_bound(triaRules_.cbegin(), triaRules_.cend(), doe, compDoe<QuadRule2D>);
  return iter != triaRules_.cend() ? *iter : *(--iter);
}

void QuadRuleManager::setTetraRule(const QuadRule3D& rule)
{
  auto insertion = tetraRules_.emplace(rule);
  if(insertion.second == false)
  {
    auto deletion = tetraRules_.erase(insertion.first);
    tetraRules_.emplace_hint(deletion, rule);
  }
}

void QuadRuleManager::setTetraRule(QuadRule3D&& rule)
{
  auto insertion = tetraRules_.emplace(std::move(rule));
  if(insertion.second == false)
  {
    auto deletion = tetraRules_.erase(insertion.first);
    tetraRules_.emplace_hint(deletion, std::move(rule));
  }
}

void QuadRuleManager::setTriaRule(const QuadRule2D& rule)
{
  auto insertion = triaRules_.emplace(rule);
  if(insertion.second == false)
  {
    auto deletion = triaRules_.erase(insertion.first);
    triaRules_.emplace_hint(deletion, rule);
  }
}

void QuadRuleManager::setTriaRule(QuadRule2D&& rule)
{
  auto insertion = triaRules_.emplace(std::move(rule));
  if(insertion.second == false)
  {
    auto deletion = triaRules_.erase(insertion.first);
    triaRules_.emplace_hint(deletion, std::move(rule));
  }
}

QuadRuleManager::QuadRuleManager()
  : faceMaps_
{
  // Here the maps from the standard 2d-simplex to the four faces of the standard
  // 3d-simplex are defined.

  // Face z = 0
  (Eigen::Matrix3d() << 0.0, 1.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0).finished(),
  // Face y = 0
  (Eigen::Matrix3d() << 1.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0).finished(),
  // Face x = 0
  (Eigen::Matrix3d() << 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        1.0, 0.0, 0.0).finished(),
  // Face x + y + z = 1
  (Eigen::Matrix3d() << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                       -1.0, -1.0, 1.0).finished()}
{
  // Here the quadrature rules over the standard 3d-simplex are initalized.
  // There are 8 rules with degree of exactness up to 8.

  // 1 point, degree of exactness = 1, [Quarteroni]
  QuadRule3D tetra1(1,{{0.25, 0.25, 0.25}}, {1./6.});

  // 4 points, degree of exactness = 2, [Quarteroni]
  QuadRule3D tetra2(2,
  {{0.585410196624968, 0.138196601125011, 0.138196601125011},
   {0.138196601125011, 0.585410196624968, 0.138196601125011},
   {0.138196601125011, 0.138196601125011, 0.585410196624968},
   {0.138196601125011, 0.138196601125011, 0.138196601125011}},
  {1./24., 1./24., 1./24., 1./24.});

  // 5 points, degree of exactness = 3, [Quarteroni]
  QuadRule3D tetra3(3,
  {{0.25, 0.25, 0.25},
   {0.5,   1./6., 1./6.},
   {1./6.,  0.5,  1./6.},
   {1./6., 1./6.,  0.5},
   {1./6., 1./6., 1./6.}},
  {-4./30., 0.075, 0.075, 0.075, 0.075});

  // 11 points, degree of exactness = 4, [Keast]
  // N.B. negative weight!
  QuadRule3D tetra4(4,
  {{0.25, 0.25, 0.25},
   {0.0714285714285714285, 0.0714285714285714285, 0.0714285714285714285},
   {0.785714285714285714, 0.0714285714285714285, 0.0714285714285714285},
   {0.0714285714285714285, 0.785714285714285714, 0.0714285714285714285},
   {0.0714285714285714285, 0.0714285714285714285, 0.785714285714285714},
   {0.399403576166799219, 0.100596423833200785, 0.100596423833200785},
   {0.100596423833200785, 0.399403576166799219, 0.100596423833200785},
   {0.100596423833200785, 0.100596423833200785, 0.399403576166799219},
   {0.100596423833200785, 0.399403576166799219, 0.399403576166799219},
   {0.399403576166799219, 0.100596423833200785, 0.399403576166799219},
   {0.399403576166799219, 0.399403576166799219, 0.100596423833200785}},
  {-0.0131555555555555550, 0.00762222222222222222, 0.00762222222222222222,
    0.00762222222222222222, 0.00762222222222222222, 0.0248888888888888880,
    0.0248888888888888880, 0.0248888888888888880,  0.0248888888888888880,
    0.0248888888888888880, 0.0248888888888888880});

  // 15 points, degree of exactness = 5, [Keast]
  QuadRule3D tetra5(5,
  {{0.25, 0.25, 0.25},
   {1./3., 1./3., 1./3.},
   {0.0, 1./3., 1./3.},
   {1./3., 0.0, 1./3.},
   {1./3., 1./3., 0.0},
   {1./11., 1./11., 1./11.},
   {8./11., 1./11., 1./11.},
   {1./11., 8./11., 1./11.},
   {1./11., 1./11., 8./11.},
   {0.4334498464263357, 0.0665501535736643, 0.0665501535736643},
   {0.0665501535736643, 0.4334498464263357, 0.0665501535736643},
   {0.0665501535736643, 0.0665501535736643, 0.4334498464263357},
   {0.0665501535736643, 0.4334498464263357, 0.4334498464263357},
   {0.4334498464263357, 0.0665501535736643, 0.4334498464263357},
   {0.4334498464263357, 0.4334498464263357, 0.0665501535736643}},
  {0.0302836780970891856,  0.00602678571428571597, 0.00602678571428571597,
   0.00602678571428571597, 0.00602678571428571597, 0.0116452490860289742,
   0.0116452490860289742,  0.0116452490860289742,  0.0116452490860289742,
   0.0109491415613864534,  0.0109491415613864534,  0.0109491415613864534,
   0.0109491415613864534,  0.0109491415613864534,  0.0109491415613864534});

  // 24 points, degree of exactness = 6, [Keast]
  QuadRule3D tetra6(6,
  {{0.214602871259151684, 0.214602871259151684, 0.214602871259151684},
   {0.356191386222544953, 0.214602871259151684, 0.214602871259151684},
   {0.214602871259151684, 0.356191386222544953, 0.214602871259151684},
   {0.214602871259151684, 0.214602871259151684, 0.356191386222544953},
   {0.0406739585346113397, 0.0406739585346113397, 0.0406739585346113397},
   {0.877978124396165982, 0.0406739585346113397, 0.0406739585346113397},
   {0.0406739585346113397, 0.877978124396165982, 0.0406739585346113397},
   {0.0406739585346113397, 0.0406739585346113397, 0.877978124396165982},
   {0.322337890142275646, 0.322337890142275646, 0.322337890142275646},
   {0.0329863295731730594, 0.322337890142275646, 0.322337890142275646},
   {0.322337890142275646, 0.0329863295731730594, 0.322337890142275646},
   {0.322337890142275646, 0.322337890142275646, 0.0329863295731730594},
   {0.0636610018750175299, 0.0636610018750175299, 0.269672331458315867},
   {0.0636610018750175299, 0.0636610018750175299, 0.603005664791649076},
   {0.0636610018750175299, 0.269672331458315867, 0.0636610018750175299},
   {0.0636610018750175299, 0.269672331458315867, 0.603005664791649076},
   {0.0636610018750175299, 0.603005664791649076, 0.0636610018750175299},
   {0.0636610018750175299, 0.603005664791649076, 0.269672331458315867},
   {0.269672331458315867, 0.0636610018750175299, 0.0636610018750175299},
   {0.269672331458315867, 0.0636610018750175299, 0.603005664791649076},
   {0.269672331458315867, 0.603005664791649076, 0.0636610018750175299},
   {0.603005664791649076, 0.0636610018750175299, 0.0636610018750175299},
   {0.603005664791649076, 0.0636610018750175299, 0.269672331458315867},
   {0.603005664791649076, 0.269672331458315867, 0.0636610018750175299}},
  {0.00665379170969464506, 0.00665379170969464506, 0.00665379170969464506,
    0.00665379170969464506, 0.00167953517588677620, 0.00167953517588677620,
    0.00167953517588677620, 0.00167953517588677620, 0.00922619692394239843,
    0.00922619692394239843, 0.00922619692394239843, 0.00922619692394239843,
    0.00803571428571428248, 0.00803571428571428248, 0.00803571428571428248,
    0.00803571428571428248, 0.00803571428571428248, 0.00803571428571428248,
    0.00803571428571428248, 0.00803571428571428248, 0.00803571428571428248,
    0.00803571428571428248, 0.00803571428571428248, 0.00803571428571428248});

  // 31 points, degree od exactness = 7, [Keast]
  // N.B. negative weight!
  QuadRule3D tetra7(7,
  {{0.5, 0.0, 0.0},
   {0.0, 0.5, 0.0},
   {0.0, 0.0, 0.5},
   {0.5, 0.5, 0.0},
   {0.5, 0.0, 0.5},
   {0.0, 0.5, 0.5},
   {0.25, 0.25, 0.25},
   {0.0782131923303186549, 0.0782131923303186549, 0.0782131923303186549},
   {0.765360423009044044, 0.0782131923303186549, 0.0782131923303186549},
   {0.0782131923303186549, 0.765360423009044044, 0.0782131923303186549},
   {0.0782131923303186549, 0.0782131923303186549, 0.765360423009044044},
   {0.121843216663904411, 0.121843216663904411, 0.121843216663904411},
   {0.634470350008286765, 0.121843216663904411, 0.121843216663904411},
   {0.121843216663904411, 0.634470350008286765, 0.121843216663904411},
   {0.121843216663904411, 0.121843216663904411, 0.634470350008286765},
   {0.332539164446420554, 0.332539164446420554, 0.332539164446420554},
   {0.00238250666073834549, 0.332539164446420554, 0.332539164446420554},
   {0.332539164446420554, 0.00238250666073834549, 0.332539164446420554},
   {0.332539164446420554, 0.332539164446420554, 0.00238250666073834549},
   {0.1, 0.1, 0.2},
   {0.1, 0.1, 0.6},
   {0.1, 0.2, 0.1},
   {0.1, 0.2, 0.6},
   {0.1, 0.6, 0.1},
   {0.1, 0.6, 0.2},
   {0.2, 0.1, 0.1},
   {0.2, 0.1, 0.6},
   {0.2, 0.6, 0.1},
   {0.6, 0.1, 0.1},
   {0.6, 0.1, 0.2},
   {0.6, 0.2, 0.1}},
  {0.000970017636684296702, 0.000970017636684296702, 0.000970017636684296702,
   0.000970017636684296702, 0.000970017636684296702, 0.000970017636684296702,
   0.0182642234661087939, 0.0105999415244141609, 0.0105999415244141609,
   0.0105999415244141609, 0.0105999415244141609, -0.0625177401143299494,
  -0.0625177401143299494, -0.0625177401143299494, -0.0625177401143299494,
   0.00489142526307353653, 0.00489142526307353653, 0.00489142526307353653,
   0.00489142526307353653, 0.0275573192239850917, 0.0275573192239850917,
   0.0275573192239850917, 0.0275573192239850917, 0.0275573192239850917,
   0.0275573192239850917, 0.0275573192239850917, 0.0275573192239850917,
   0.0275573192239850917, 0.0275573192239850917, 0.0275573192239850917,
   0.0275573192239850917});

  // 45 points, degree of exactness = 8, [Keast]
  // N.B. negative weight!
  QuadRule3D tetra8(8,
  {{0.25, 0.25, 0.25},
   {0.127470936566639015, 0.127470936566639015, 0.127470936566639015},
   {0.617587190300082967, 0.127470936566639015, 0.127470936566639015},
   {0.127470936566639015, 0.617587190300082967, 0.127470936566639015},
   {0.127470936566639015, 0.127470936566639015, 0.617587190300082967},
   {0.0320788303926322960, 0.0320788303926322960, 0.0320788303926322960},
   {0.903763508822103123, 0.0320788303926322960, 0.0320788303926322960},
   {0.0320788303926322960, 0.903763508822103123, 0.0320788303926322960},
   {0.0320788303926322960, 0.0320788303926322960, 0.903763508822103123},
   {0.450222904356718978, 0.0497770956432810185, 0.0497770956432810185},
   {0.0497770956432810185, 0.450222904356718978, 0.0497770956432810185},
   {0.0497770956432810185, 0.0497770956432810185, 0.450222904356718978},
   {0.0497770956432810185, 0.450222904356718978, 0.450222904356718978},
   {0.450222904356718978, 0.0497770956432810185, 0.450222904356718978},
   {0.450222904356718978, 0.450222904356718978, 0.0497770956432810185},
   {0.316269552601450060, 0.183730447398549945, 0.183730447398549945},
   {0.183730447398549945, 0.316269552601450060, 0.183730447398549945},
   {0.183730447398549945, 0.183730447398549945, 0.316269552601450060},
   {0.183730447398549945, 0.316269552601450060, 0.316269552601450060},
   {0.316269552601450060, 0.183730447398549945, 0.316269552601450060},
   {0.316269552601450060, 0.316269552601450060, 0.183730447398549945},
   {0.0229177878448171174, 0.231901089397150906, 0.231901089397150906},
   {0.231901089397150906, 0.0229177878448171174, 0.231901089397150906},
   {0.231901089397150906, 0.231901089397150906, 0.0229177878448171174},
   {0.513280033360881072, 0.231901089397150906, 0.231901089397150906},
   {0.231901089397150906, 0.513280033360881072, 0.231901089397150906},
   {0.231901089397150906, 0.231901089397150906, 0.513280033360881072},
   {0.0229177878448171174, 0.513280033360881072, 0.231901089397150906},
   {0.0229177878448171174, 0.231901089397150906, 0.513280033360881072},
   {0.513280033360881072, 0.0229177878448171174, 0.231901089397150906},
   {0.513280033360881072, 0.231901089397150906, 0.0229177878448171174},
   {0.231901089397150906, 0.0229177878448171174, 0.513280033360881072},
   {0.231901089397150906, 0.513280033360881072, 0.0229177878448171174},
   {0.730313427807538396, 0.0379700484718286102, 0.0379700484718286102},
   {0.0379700484718286102, 0.730313427807538396, 0.0379700484718286102},
   {0.0379700484718286102, 0.0379700484718286102, 0.730313427807538396},
   {0.193746475248804382, 0.0379700484718286102, 0.0379700484718286102},
   {0.0379700484718286102, 0.193746475248804382, 0.0379700484718286102},
   {0.0379700484718286102, 0.0379700484718286102, 0.193746475248804382},
   {0.0379700484718286102, 0.730313427807538396, 0.193746475248804382},
   {0.0379700484718286102, 0.193746475248804382, 0.730313427807538396},
   {0.730313427807538396, 0.193746475248804382, 0.0379700484718286102},
   {0.730313427807538396, 0.0379700484718286102, 0.193746475248804382},
   {0.193746475248804382, 0.730313427807538396, 0.0379700484718286102},
   {0.193746475248804382, 0.0379700484718286102, 0.730313427807538396}},
  {-0.0393270066412926145, 0.00408131605934270525, 0.00408131605934270525,
   0.00408131605934270525, 0.00408131605934270525, 0.000658086773304341943,
   0.000658086773304341943, 0.000658086773304341943, 0.000658086773304341943,
   0.00438425882512284693, 0.00438425882512284693, 0.00438425882512284693,
   0.00438425882512284693, 0.00438425882512284693, 0.00438425882512284693,
   0.0138300638425098166, 0.0138300638425098166, 0.0138300638425098166,
   0.0138300638425098166, 0.0138300638425098166, 0.0138300638425098166,
   0.00424043742468372453, 0.00424043742468372453, 0.00424043742468372453,
   0.00424043742468372453, 0.00424043742468372453, 0.00424043742468372453,
   0.00424043742468372453, 0.00424043742468372453, 0.00424043742468372453,
   0.00424043742468372453, 0.00424043742468372453, 0.00424043742468372453,
   0.00223873973961420164, 0.00223873973961420164, 0.00223873973961420164,
   0.00223873973961420164, 0.00223873973961420164, 0.00223873973961420164,
   0.00223873973961420164, 0.00223873973961420164, 0.00223873973961420164,
   0.00223873973961420164, 0.00223873973961420164, 0.00223873973961420164});

    // Here the quadrature rules over the standard 2d-simplex are defined.
    // There are 10 rules with degree of exactness up to 10.

    // 1 point, degree of exactness = 1, [Dunavant]
    QuadRule2D tria1(1,
    {{1./3., 1./3.}},
    {0.5});

    // 3 points, degree of exactness = 2, [Dunavant]
    QuadRule2D tria2(2,
    {{2./3., 1./6.},
     {1./6., 2./3.},
     {1./6., 1./6.}},
    {1./6., 1./6., 1./6.});

    // 4 points, degree of exactness = 3, [Dunavant]
    // N.B. negative weight!
    QuadRule2D tria3(3,
    {{1./3., 1./3.},
     {0.2, 0.2},
     {0.6, 0.2},
     {0.2, 0.6}},
    {-9./32., 25./96., 25./96., 25./96.});

    // 6 points, degree of exactness = 4, [Dunavant]
    QuadRule2D tria4(4,
    {{0.108103018168070, 0.445948490915965},
     {0.445948490915965, 0.445948490915965},
     {0.445948490915965, 0.108103018168070},
     {0.816847572980459, 0.091576213509771},
     {0.091576213509771, 0.091576213509771},
     {0.091576213509771, 0.816847572980459}},
    {0.1116907948390055, 0.1116907948390055, 0.1116907948390055,
     0.054975871827661, 0.054975871827661, 0.054975871827661});

    // 7 points, degree of exactness = 5, [Dunavant]
    QuadRule2D tria5(5,
    {{1./3., 1./3.},
     {0.470142064105115, 0.470142064105115},
     {0.059715871789770, 0.470142064105115},
     {0.470142064105115, 0.059715871789770},
     {0.101286507323456, 0.101286507323456},
     {0.797426985353087, 0.101286507323456},
     {0.101286507323456, 0.797426985353087}},
    {9./80., 0.066197076394253, 0.066197076394253, 0.066197076394253, 0.0629695902724135,
     0.0629695902724135, 0.0629695902724135});

    // 12 points, degree of exactness = 6, [Dunavant]
    QuadRule2D tria6(6,
    {{0.501426509658179, 0.249286745170910},
     {0.249286745170910, 0.249286745170910},
     {0.249286745170910, 0.501426509658179},
     {0.873821971016996, 0.063089014491502},
     {0.063089014491502, 0.063089014491502},
     {0.063089014491502, 0.873821971016996},
     {0.053145049844817, 0.310352451033784},
     {0.310352451033784, 0.053145049844817},
     {0.053145049844817, 0.636502499121399},
     {0.636502499121399, 0.053145049844817},
     {0.636502499121399, 0.310352451033784},
     {0.310352451033784, 0.636502499121399}},
    {0.0583931378631895, 0.0583931378631895, 0.0583931378631895, 0.0254224531851035,
     0.0254224531851035, 0.0254224531851035, 0.041425537809187, 0.041425537809187,
     0.041425537809187, 0.041425537809187, 0.041425537809187, 0.041425537809187});

    // 13 points, degree of exactness = 7, [Dunavant]
    // N.B. negative weight!
    QuadRule2D tria7(7,
    {{1./3., 1./3.},
     {0.479308067841920, 0.260345966079040},
     {0.260345966079040, 0.260345966079040},
     {0.260345966079040, 0.479308067841920},
     {0.869739794195568, 0.065130102902216},
     {0.065130102902216, 0.065130102902216},
     {0.065130102902216, 0.869739794195568},
     {0.048690315425316, 0.312865496004876},
     {0.312865496004876, 0.048690315425316},
     {0.048690315425316, 0.638444188569810},
     {0.638444188569810, 0.048690315425316},
     {0.638444188569810, 0.312865496004874},
     {0.312865496004874, 0.638444188569810}},
    {-0.074785022233841, 0.087807628716604, 0.087807628716604, 0.087807628716604,
      0.026673617804419, 0.026673617804419, 0.026673617804419, 0.0385568804451285,
      0.0385568804451285, 0.0385568804451285, 0.0385568804451285, 0.0385568804451285,
      0.0385568804451285});

    // 16 points, degree of exactness = 8, [Dunavant]
    QuadRule2D tria8(8,
    {{1./3., 1./3.},
     {0.459292588292723, 0.459292588292723},
     {0.459292588292723, 0.081414823414554},
     {0.081414823414554, 0.459292588292723},
     {0.170569307751760, 0.170569307751760},
     {0.170569307751760, 0.658861384496480},
     {0.658861384496480, 0.170569307751760},
     {0.050547228317031, 0.050547228317031},
     {0.050547228317031, 0.898905543365938},
     {0.898905543365938, 0.050547228317031},
     {0.263112829634638, 0.728492392955404},
     {0.728492392955404, 0.008394777409958},
     {0.008394777409958, 0.263112829634638},
     {0.728492392955404, 0.263112829634638},
     {0.263112829634638, 0.008394777409958},
     {0.008394777409958, 0.728492392955404}},
    {0.0721578038388935, 0.0475458171336425, 0.0475458171336425, 0.0475458171336425,
     0.051608685267359, 0.051608685267359, 0.051608685267359, 0.016229248811599,
     0.016229248811599, 0.016229248811599, 0.0136151570872175, 0.0136151570872175,
     0.0136151570872175, 0.0136151570872175, 0.0136151570872175, 0.0136151570872175});

    // 19 points, degree of exactness = 9, [Dunavant]
    QuadRule2D tria9(9,
    {{1./3., 1./3.},
     {0.020634961602525, 0.489682519198738},
     {0.489682519198738, 0.020634961602525},
     {0.489682519198738, 0.489682519198738},
     {0.125820817014127, 0.437089591492937},
     {0.437089591492937, 0.125820817014127},
     {0.437089591492937, 0.437089591492937},
     {0.623592928761935, 0.188203535619033},
     {0.188203535619033, 0.623592928761935},
     {0.188203535619033, 0.188203535619033},
     {0.910540973211095, 0.044729513394453},
     {0.044729513394453, 0.910540973211095},
     {0.044729513394453, 0.044729513394453},
     {0.036838412054736, 0.221962989160766},
     {0.036838412054736, 0.741198598784498},
     {0.221962989160766, 0.036838412054736},
     {0.221962989160766, 0.741198598784498},
     {0.741198598784498, 0.036838412054736},
     {0.741198598784498, 0.221962989160766}},
    {0.0485678981413995, 0.0156673501135695, 0.0156673501135695, 0.0156673501135695,
     0.038913770502387, 0.038913770502387, 0.038913770502387, 0.039823869463605,
     0.039823869463605, 0.039823869463605, 0.012788837829349, 0.012788837829349,
     0.012788837829349, 0.0216417696886445, 0.0216417696886445, 0.0216417696886445,
     0.0216417696886445, 0.0216417696886445, 0.0216417696886445});

    // 25 points, degree of exactness = 10, [Dunavant]
    QuadRule2D tria10(10,
    {{1./3., 1./3.},
     {0.028844733232685, 0.485577633383657},
     {0.485577633383657, 0.028844733232685},
     {0.485577633383657, 0.485577633383657},
     {0.781036849029926, 0.109481575485037},
     {0.109481575485037, 0.781036849029926},
     {0.109481575485037, 0.109481575485037},
     {0.141707219414880, 0.307939838764121},
     {0.141707219414880, 0.550352941820999},
     {0.307939838764121, 0.550352941820999},
     {0.307939838764121, 0.141707219414880},
     {0.550352941820999, 0.307939838764121},
     {0.550352941820999, 0.141707219414880},
     {0.025003534762686, 0.246672560639903},
     {0.025003534762686, 0.728323904597411},
     {0.246672560639903, 0.728323904597411},
     {0.246672560639903, 0.025003534762686},
     {0.728323904597411, 0.246672560639903},
     {0.728323904597411, 0.025003534762686},
     {0.009540815400299, 0.066803251012200},
     {0.009540815400299, 0.923655933587500},
     {0.066803251012200, 0.009540815400299},
     {0.066803251012200, 0.923655933587500},
     {0.923655933587500, 0.066803251012200},
     {0.923655933587500, 0.009540815400299}},
    {0.045408995191377, 0.0183629788782335, 0.0183629788782335, 0.0183629788782335,
     0.022660529717764, 0.022660529717764, 0.022660529717764, 0.036378958422710,
     0.036378958422710, 0.036378958422710, 0.036378958422710, 0.036378958422710,
     0.036378958422710, 0.0141636212655285, 0.0141636212655285, 0.0141636212655285,
     0.0141636212655285, 0.0141636212655285,0.0141636212655285, 0.0047108334818665,
     0.0047108334818665, 0.0047108334818665, 0.0047108334818665, 0.0047108334818665,
     0.0047108334818665});

  setTetraRule(std::move(tetra1));
  setTetraRule(std::move(tetra2));
  setTetraRule(std::move(tetra3));
  setTetraRule(std::move(tetra4));
  setTetraRule(std::move(tetra5));
  setTetraRule(std::move(tetra6));
  setTetraRule(std::move(tetra7));
  setTetraRule(std::move(tetra8));

  setTriaRule(std::move(tria1));
  setTriaRule(std::move(tria2));
  setTriaRule(std::move(tria3));
  setTriaRule(std::move(tria4));
  setTriaRule(std::move(tria5));
  setTriaRule(std::move(tria6));
  setTriaRule(std::move(tria7));
  setTriaRule(std::move(tria8));
  setTriaRule(std::move(tria9));
  setTriaRule(std::move(tria10));
}

} // namespace PolyDG
