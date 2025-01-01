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
  using B       = Soyuz21b_Stage1_Booster<'B'>;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  // Thrusts and Times:
  cout << "# BoosterMaxFullThrDur  : " << B::MaxFullThrustDur         << endl;
  cout << "# BoosterCutOffTime     : " << Time(B::CutOffTime)         << endl;
  cout << "# BoosterMaxBurnTime    : " << Time(B::MaxBurnTime)        << endl;
  cout << "# BoosterThrustVac      : " << B::ThrustEngVac   / g0      << endl;
  cout << "# BoosterThrustVernSL1  : " << B::ThrustVernSL1  / g0      << endl;
  cout << "# BoosterThrustVernVac1 : " << B::ThrustVernVac1 / g0      << endl;
  cout << "# BoosterThrustMainSL   : " << B::ThrustMainSL   / g0      << endl;
  cout << "# BoosterThrustMainVac  : " << S2::ThrustMainVac / g0      << endl;

  // Max and Real Oxid and Fuel Loads:
  cout << "# BoosterMaxOxidLoad    : " << B::OxidTankMC               << endl;
  cout << "# BoosterActOxidLoad    : " << B::OxidLoadRatio  * 100.0   << " %"
       << endl;
  cout << "# BoosterMaxFuelLoad    : " << B::FuelTankMC               << endl;
  cout << "# BoosterActFuelLoad    : " << B::FuelLoadRatio  * 100.0   << " %"
       << endl;
  cout << "# BoosterActLiqN2Load   : " << B::LiqN2LoadRatio * 100.0   << " %"
       << endl;
  cout << "# BoosterActH2O2Load    : " << B::H2O2LoadRatio  * 100.0   << " %"
       << endl;

  // Geometry:
  cout << "# BoosterLen            : " << B::H                        << endl;
  cout << "# BoosterAlpha          : "
       << To_Angle_deg(Angle(ATan(B ::TanAlpha)))                     << endl;
  cout << "# BoosterAlphaTop       : "
       << To_Angle_deg(Angle(ATan(B ::TanAlphaTop)))                  << endl;
  cout << "# BoosterAlphaMid       : "
       << To_Angle_deg(Angle(ATan(S2::TanAlphaMid)))                  << endl;
  cout << "# BoosterOxidTankMidD   : " << B::OxidTankMidD             << endl;
  cout << "# BoosterInterTankLoD   : " << B::InterTankLoD             << endl;
  cout << "# BoosterOxidFuelTankGap: " << B::OxidFuelTankGap          << endl;
  cout << "# BoosterFuelTankLowX   : " << B::FuelTankBtm .GetLow()[0] << endl;
  cout << "# BoosterLiqN2TankUpX   : " << B::LiqN2TankTop.GetUp ()[0] << endl;
  cout << "# BoosterLiqN2TankLowX  : " << B::LiqN2TankBtm.GetLow()[0] << endl;
  cout << '#' << endl;
  cout << "# BoosterH2O2TankUpX    : " << B::H2O2TankTop .GetUp ()[0] << endl;
  cout << "# BoosterH2O2TankLowX   : " << B::H2O2TankBtm .GetLow()[0] << endl;
  cout << "# BoosterTailCylLowX    : " << B::TailLowX                 << endl;
  cout << "# BoosterTailGapH       : " << B::TailGapH                 << endl;
  cout << "# BoosterNozzlesExtH    : " << B::NozzlesExtH              << endl;
  cout << "# BoosterNozzlesLowX    : " << B::NozzlesLowX              << endl;

  // Header for the following table:
  cout << '#' << endl;
  cout << "# (1)Time\t(2)TotalMass\t(3)FuelMass\t(4)OxidMass\t(5)CoM_x\t"
          "(6)J_x\t(7)J_y\t(8)J_z\t(9)JDot_x\t(10)JDot_y\t(11)JDot_z" << endl;
  cout << '#' << endl;
  return 0;
}
