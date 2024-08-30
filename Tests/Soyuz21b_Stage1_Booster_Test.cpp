// vim:ts=2:et
//===========================================================================//
//                "Tests/Soyuz21b_Stage1_Booster_Test.cpp":                  //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#include "SpaceBallistics/Utils.hpp"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  using      B  = Soyuz21b_Stage1_Booster<'B'>;
  using      S2 = Soyuz21b_Stage2;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  // Thrusts and Times:
  cout << "# BoosterMaxFullThrDur  : " << B::MaxFullThrustDur       << endl;
  cout << "# BoosterCutOffTime     : " << B::CutOffTime             << endl;
  cout << "# BoosterMaxBurnTime    : " << B::MaxBurnTime            << endl;
  cout << "# BoosterThrustVac      : " << B::ThrustEngVac   / g0    << endl;
  cout << "# BoosterThrustVernSL1  : " << B::ThrustVernSL1  / g0    << endl;
  cout << "# BoosterThrustVernVac1 : " << B::ThrustVernVac1 / g0    << endl;
  cout << "# BoosterThrustMainSL   : " << B::ThrustMainSL   / g0    << endl;
  cout << "# BoosterThrustMainVac  : " << S2::ThrustMainVac / g0    << endl;

  // Max and Real Oxid and Fuel Loads:
  cout << "# BoosterMaxOxidLoad    : " << B::OxidTankMC             << endl;
  cout << "# BoosterActOxidLoad    : " << B::OxidLoadRatio * 100.0  << " %"
       << endl;
  cout << "# BoosterMaxFuelLoad    : " << B::FuelTankMC             << endl;
  cout << "# BoosterActFuelLoad    : " << B::FuelLoadRatio * 100.0  << " %"
       << endl;

  // Geometry:
  cout << "# BoosterLen            : " << B::H                      << endl;
  cout << "# BoosterAlpha          : "
       << To_Angle_deg(Angle(ATan(B ::TanAlpha)))                   << endl;
  cout << "# BoosterAlphaTop       : "
       << To_Angle_deg(Angle(ATan(B ::TanAlphaTop)))                << endl;
  cout << "# BoosterAlphaMid       : "
       << To_Angle_deg(Angle(ATan(S2::TanAlphaMid)))                << endl;
  cout << "# BoosterOxidTankMidD   : " << B::OxidTankMidD           << endl;
  cout << "# BoosterInterTankLoD   : " << B::InterTankLoD           << endl;
  cout << "# BoosterOxidFuelTankGap: " << B::OxidFuelTankGap        << endl;
  cout << "# BoosterTailCylLowX    : " << B::TailLowX               << endl;

  // Header for the following table:
  cout << '#' << endl;
  cout << "# Time\tTotalMass\tFuelMass\tOxidMass\tCoM_x\tJ_x\tJ_y" << endl;
  cout << '#' << endl;
  return 0;
}
