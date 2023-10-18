// vim:ts=2:et
//===========================================================================//
//                            "Tests/Test0.cpp":                             //
//===========================================================================//
#include "SpaceBallistics/Soyuz21b.h"
#include "SpaceBallistics/Soyuz21b.hpp"
#include <iostream>

int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  using namespace Soyuz21b_Consts;

  // The LV Object:
  Soyuz21b lv;

  // Gap in Stage3 between the bottom of FuelTank and the top of OxidTank:
  Len  fuelTankBot3 = Stage3FuelTankLoX0 + Stage3FuelTankLoH;
  Len  oxidTankTop3 = Stage3OxidTankUpX0 - Stage3OxidTankUpH;
  cout << "Stage3FuelOxidGap   : " << ToStr(oxidTankTop3 - fuelTankBot3).data()
       << endl;

  // The max possible fuel and oxidiser load for Stage3 (using the max geometr-
  // ical tank volumes):
  Mass maxFuelMass3 = Stage3FuelTankVol * RG1Dens;
  Mass maxOxidMass3 = Stage3OxidTankVol * LOXDens;
  Time maxBurnTime3 = Stage3SpentMass / Stage3MassRate;

  cout << "Stage3MaxFuelLoad   : " << ToStr(maxFuelMass3).data()   << endl;
  cout << "Stage3ActFuelLoad   : "
       << (double(Stage3FuelMass / maxFuelMass3) * 100.0) << " %"  << endl;

  cout << "Stage3MaxOxidLoad   : " << ToStr(maxOxidMass3).data()   << endl;
  cout << "Stage3ActOxidLoad   : "
       << (double(Stage3OxidMass / maxOxidMass3) * 100.0) << " %"  << endl;

  cout << "Stage3MassRate      : " << ToStr(Stage3MassRate).data() << endl;
  cout << "Stage3MaxBurnTime   : " << ToStr(maxBurnTime3).data()   << endl;

  // Percentage of the StaticThrust in the over-all thrust of Stage3:
  cout << "Stage3StaticThrust  : " << ToStr(Stage3StaticThrust).data()
       << "\n\t\t\t(" << (double(Stage3StaticThrust/Stage3ThrustVac) * 100.0)
       << " % of over-all)" << endl;

  // Length and Mass of the Stage3Aft:
  cout << "Stage3AftLen        : " << ToStr(Stage3AftH).data()  << endl;
  cout << "Stage3AftMass       : "
       << ToStr(lv.m_stage3EmptyMass - lv.m_stage3NoAftEmptyMass).data()
       << endl;
  // MoI:
  cout << "Stage3EmptyMoIY     : " << ToStr(lv.m_stage3EmptyMoIY).data()
       << endl;
  cout << "Stage3NoAftEmptyMoIY: " << ToStr(lv.m_stage3NoAftEmptyMoIY).data()
       << endl;

  return 0;
}
