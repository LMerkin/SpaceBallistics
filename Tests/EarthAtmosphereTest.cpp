// vim:ts=2:et
//===========================================================================//
//                     "Tests/EarthAtmosphereTest.cpp":                      //
//===========================================================================//
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include <iostream>

int main()
{
  namespace EAM = SpaceBallistics::EarthAtmosphereModel;
  using namespace SpaceBallistics;
  using namespace std;

  for (LenK z = 0.0_km; z <= 1201.0_km; z += 1.0_km)
  {
    auto [p, rho, T, a] = EAM::GetAtmConds(z);
    cout << z.Magnitude() << '\t' << T  .Magnitude() << '\t'
         << p.Magnitude() << '\t' << rho.Magnitude() << '\t'
         << a.Magnitude() << std::endl;
  }
  return 0;
}
