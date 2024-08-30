// vim:ts=2:et
//===========================================================================//
//                     "Tests/Soyuz21b_Stage2_Test.cpp":                     //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  using      S2 = Soyuz21b_Stage2;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  // Thrusts and Times:
  cout << "# Stage2MaxFullThrDur : " << S2::MaxFullThrustDur       << endl;
  cout << "# Stage2CutOffTime    : " << S2::CutOffTime             << endl;
  cout << "# Stage2MaxBurnTime   : " << S2::MaxBurnTime            << endl;
  cout << "# Stage2ThrustSL      : " << S2::ThrustEngSL    / g0    << endl;
  cout << "# Stage2ThrustVac     : " << S2::ThrustEngVac   / g0    << endl;
  cout << "# Stage2ThrustVernSL1 : " << S2::ThrustVernSL1  / g0    << endl;
  cout << "# Stage2ThrustVernVac1: " << S2::ThrustVernVac1 / g0    << endl;
  cout << "# Stage2ThrustMainSL  : " << S2::ThrustMainSL   / g0    << endl;
  cout << "# Stage2ThrustMainVac : " << S2::ThrustMainVac  / g0    << endl;

  // Max and Real Oxid and Fuel Loads:
  cout << "# Stage2MaxOxidLoad   : " << S2::OxidTankMC             << endl;
  cout << "# Stage2ActOxidLoad   : " << S2::OxidLoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxFuelLoad   : " << S2::FuelTankMC             << endl;
  cout << "# Stage2ActFuelLoad   : " << S2::FuelLoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxH2O2Load   : " << S2::H2O2TankMC             << endl;
  cout << "# Stage2ActH2O2Load   : " << S2::H2O2LoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxN2Load     : " << S2::N2TankMC               << endl;
  cout << "# Stage2ActN2Load     : " << S2::N2LoadRatio   * 100.0  << " %"
       << endl;

  // Remnants:
  cout << "# Stage2Fuel@CutOff   : " << S2::FuelMassC << endl;
  cout << "# Stage2Oxid@CutOff   : " << S2::OxidMassC << endl;

  // Geometry:
  cout << "# Stage2OverAllLen    : " << S2::H                      << endl;
  cout << "# Stage2EngineCoMX    : " << S2::EngineCoMX             << endl;
  cout << "# Stage2TailEnclLowX  : " << S2::TailEnclLowX           << endl;
  cout << "# Stage2NozzlesLowX   : " << S2::NozzlesLowX            << endl;

  // Header for the following table:
  cout << '#' << endl;
  cout << "# Time\tTotalMass\tFuelMass\tOxidMass\tCoM_x\tJ_x\tJ_y" << endl;
  cout << '#' << endl;

  //-------------------------------------------------------------------------//
  // Stage2 Params as a function of Flight Time:                             //
  //-------------------------------------------------------------------------//
  // XXX: Currently assuming no Thrust Vector Control, and Sea-Level Atmospheric
  // Pressure:
  //
  S2::VernDeflections vernDefls0;      // All 0s by default 

  for (Time t = 0.0_sec; t <= 300.0_sec; t += 0.1_sec)
  {
    StageDynParams<LVSC::Soyuz21b> dp = S2::GetDynParams(t, p0, vernDefls0);

    assert(IsZero(dp.m_com[1]) && IsZero(dp.m_com[2]) &&
           dp.m_mois[1]        == dp.m_mois[2]);

    cout << t.Magnitude  ()           << '\t'
         << dp.m_fullMass.Magnitude() << '\t'
         << dp.m_fuelMass.Magnitude() << '\t'
         << dp.m_oxidMass.Magnitude() << '\t'
         << dp.m_com [0] .Magnitude() << '\t'
         << dp.m_mois[0] .Magnitude() << '\t'
         << dp.m_mois[1] .Magnitude() << endl;
  }
  cout << "# Stage2FullMR        : " << S2::FullMR     << endl;
  cout << "# Stage2MinEndMass    : " << S2::MinEndMass << endl;
  return 0;
}
