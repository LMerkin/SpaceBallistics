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
  namespace  SC = Soyuz21b_Consts;
  using      S2 = Soyuz21b_Stage2;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  cout << "# Stage2PropMargin  : "
       << (double(S2::MaxFlightTime / S2::CutOffTime) - 1.0) * 100.0
       << " %" << endl;

  // Max and Real Oxid and Fuel Loads:
  cout << "# Stage2MaxOxidLoad : " << S2::OxidTankMC             << endl;
  cout << "# Stage2ActOxidLoad : " << S2::OxidLoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxFuelLoad : " << S2::FuelTankMC             << endl;
  cout << "# Stage2ActFuelLoad : " << S2::FuelLoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxH2O2Load : " << S2::H2O2TankMC             << endl;
  cout << "# Stage2ActH2O2Load : " << S2::H2O2LoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage2MaxN2Load   : " << S2::N2TankMC               << endl;
  cout << "# Stage2ActN2Load   : " << S2::N2LoadRatio   * 100.0  << " %"
       << endl;

  // Over-All Length:
  cout << "# Stage2OverAllLen  : " << S2::H                      << endl;
  cout << "# EngineNozzlesLow1 : " << S2::EngineNozzlesLow1      << endl;
  cout << "# EngineNozzlesLow2 : " << S2::EngineNozzlesLow2      << endl;

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

//  assert(IsZero(dp.m_com[1]) && IsZero(dp.m_com[2]) &&
//         dp.m_mois[1]        == dp.m_mois[2]);

    cout << t.Magnitude  ()           << '\t'
         << dp.m_fullMass.Magnitude() << '\t'
         << dp.m_fuelMass.Magnitude() << '\t'
         << dp.m_oxidMass.Magnitude() << '\t'
//       << dp.m_com [0] .Magnitude() << '\t'
//       << dp.m_mois[0] .Magnitude() << '\t'
//       << dp.m_mois[1] .Magnitude()
         << endl;
  }
  cout << "# FullMR            : " << S2::FullMR     << endl;
  cout << "# MinEndMass        : " << S2::MinEndMass << endl;
  return 0;
}
