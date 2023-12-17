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

  //-------------------------------------------------------------------------//
  // MoIs and CoM of a Cylinder:                                             //
  //-------------------------------------------------------------------------//
  constexpr Len     D   = 2.0_m;
  constexpr Len     h   = 5.0_m;
  constexpr Density rho(1.0);
  constexpr Len     R   = D / 2.0;
  constexpr Vol     v   = Pi<double> * Sqr(R) * h;
  constexpr Mass    propMass = v * rho;

  constexpr TrCone  cyl {0.0_m, D, D, h, rho, 1.0_kg};
  Len com [3];
  MoI mois[3];
  cyl(propMass, false, com, mois);

  // Diffs:
  Len dC  [3] {Abs(com[0] - h/2.0), Abs(com[1]), Abs(com[2])};

  cout << "CYLINDER:" << endl;
  cout << "dxC=" << dC[0] << ", dyC=" << dC[1] << ", dzC=" << dC[2] << endl;

  // Expected MoIs for the Cylinder:
  MoI JxE  = propMass *  Sqr(R) / 2.0;
  MoI JyzE = propMass * (Sqr(R) / 4.0 + Sqr(h) / 3.0);

  cout <<   "dJx=" << Abs(mois[0] - JxE)  << ", dJy=" << Abs(mois[1] - JyzE)
       << ", dJz=" << Abs(mois[2] - JyzE) << endl;

  return 0;
}
