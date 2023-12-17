// vim:ts=2:et
//===========================================================================//
//                                "CETest.cpp":                              //
//                   Tests of the "ConstrElement" Class                      //
//===========================================================================//
#include "SpaceBallistics/ConstrElement.hpp"
#include <iostream>

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  // MoIs and CoM of a Cylinder:
  constexpr Len     D   = 2.0_m;
  constexpr Len     h   = 5.0_m;
  constexpr Density rho(1.0);
  constexpr Len     R   = D / 2.0;
  constexpr Vol     v   = Pi<double> * Sqr(R) * h;
  constexpr Mass    propMass = v * rho;

  TrCone cyl(0.0_m, D, D, h, rho, 1.0_kg);

  Len com [3];
  MoI mois[3];
  cyl(propMass, false, com, mois);

  cout << "V="    << v      << endl;
  cout << "  xC=" << com[0] << ", yC=" << com[1] << ", zC="   << com[2]
       << endl;
  cout << "Expected: xC="   << (h/2.0) << ", yC=zC=" << 0.0_m << endl << endl;

  cout <<   "Jx=" << mois[0]
       << ", Jy=" << mois[1]
       << ", Jz=" << mois[2] << endl;
  // Expected MoIs for the Cylinder:
  MoI JxE = propMass *  Sqr(R) / 2.0;
  MoI JyE = propMass * (Sqr(R) / 4.0 + Sqr(h) / 3.0);

  cout << "Expected: Jx=" << JxE << ", Jy=Jz=" << JyE << std::endl;
  return 0;
}
