// vim:ts=2:et
//===========================================================================//
//                            "Tests/Test0.cpp":                             //
//===========================================================================//
#include "Soyuz_21b.hpp"
#include <iostream>

using namespace SpaceBallistics;
using namespace Soyuz_21b_Consts;
using namespace std;

namespace
{
  struct BoilerPlate {};
}

int main()
{
  Soyuz_21b<BoilerPlate> rks { BoilerPlate() };

  // Compute the total fuel and oxidiser load for Stage3:
  Mass_kg  maxFuelMass3 = Stage3FuelTankVol * RG1Dens;
  Mass_kg  maxOxidMass3 = Stage3OxidTankVol * LOXDens;
  Time_sec burnTime3    = (Stage3FuelMass + Stage3OxidMass) / Stage3MassRate;

  cout << "Stage3MaxFuelLoad: " << ToStr(maxFuelMass3).data()   << endl;
  cout << "Stage3MaxOxidLoad: " << ToStr(maxOxidMass3).data()   << endl;
  cout << "Stage3MassRate   : " << ToStr(Stage3MassRate).data() << endl;
  cout << "Stage3BurnTime   : " << ToStr(burnTime3).data()      << endl;

  return 0;
}
