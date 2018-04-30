#include "LegendreBasis.hpp"

namespace PolyDG
{

template<>
Real LegendreBasis<0>(Real /*x*/)
{
  return 1.0;
}

template<>
Real LegendreBasis<1>(Real x)
{
  return x;
}

template<>
Real LegendreBasis<2>(Real x)
{
  return 1.5 * pow(x, 2) - 0.5;
}

template<>
Real LegendreBasis<3>(Real x)
{
  return x * (2.5 * pow(x, 2) - 1.5);
}

template<>
Real LegendreBasis<4>(Real x)
{
  return (35.0 * pow(x, 4) - 30.0 * pow(x, 2) + 3.0) / 8.0;
}

template<>
Real LegendreBasis<5>(Real x)
{
  return (63.0 * pow(x, 5) - 70.0 * pow(x, 3) + 15.0 * x) / 8.0;
}

template<>
Real LegendreBasis<6>(Real x)
{
  return (231.0 * pow(x, 6) - 315.0 * pow(x, 4) + 105.0 * pow(x, 2) - 5.0) / 16.0;
}

template<>
Real LegendreBasisDer<0>(Real /*x*/)
{
  return 0.0;
}

template<>
Real LegendreBasisDer<1>(Real /*x*/)
{
  return 1.0;
}

template<>
Real LegendreBasisDer<2>(Real x)
{
  return 3.0 * x;
}

template<>
Real LegendreBasisDer<3>(Real x)
{
  return 7.5 * pow(x, 2) - 3.0;
}

template<>
Real LegendreBasisDer<4>(Real x)
{
  return (140.0 * pow(x, 3) - 60.0 * x) / 8.0;
}

template<>
Real LegendreBasisDer<5>(Real x)
{
  return (315.0 * pow(x, 4) - 210.0 * pow(x, 2) + 15.0) / 8.0;
}

template<>
Real LegendreBasisDer<6>(Real x)
{
  return (1386.0 * pow(x, 5) - 1260.0 * pow(x, 3) + 210.0 * x) / 16.0;
}

} // namespace PolyDG
