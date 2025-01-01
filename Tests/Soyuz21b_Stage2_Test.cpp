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
  cout << "# Stage2CutOffTime    : " << Time(S2::CutOffTime)       << endl;
  cout << "# Stage2MaxBurnTime   : " << Time(S2::MaxBurnTime)      << endl;
  cout << "# Stage2ThrustSL      : " << S2::ThrustEngSL    / g0    << endl;
  cout << "# Stage2ThrustVac     : " << S2::ThrustEngVac   / g0    << endl;
  cout << "# Stage2ThrustVernSL1 : " << S2::ThrustVernSL1  / g0    << endl;
  cout << "# Stage2ThrustVernVac1: " << S2::ThrustVernVac1 / g0    << endl;
  cout << "# Stage2ThrustMainSL  : " << S2::ThrustMainSL   / g0    << endl;
  cout << "# Stage2ThrustMainVac : " << S2::ThrustMainVac  / g0    << endl;

  // Max and Real Oxid and Fuel Loads:
  cout << "# Stage2MaxOxidLoad   : " << S2::OxidTankMC             << endl;
  cout << "# Stage2ActOxidLoad   : " << S2::OxidLoadRatio  * 100.0 << " %"
       << endl;
  cout << "# Stage2MaxFuelLoad   : " << S2::FuelTankMC             << endl;
  cout << "# Stage2ActFuelLoad   : " << S2::FuelLoadRatio  * 100.0 << " %"
       << endl;
  cout << "# Stage2MaxH2O2Load   : " << S2::H2O2TankMC             << endl;
  cout << "# Stage2ActH2O2Load   : " << S2::H2O2LoadRatio  * 100.0 << " %"
       << endl;
  cout << "# Stage2MaxN2Load     : " << S2::LiqN2TankMC            << endl;
  cout << "# Stage2ActN2Load     : " << S2::LiqN2LoadRatio * 100.0 << " %"
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
  cout << "# (1)Time\t(2)TotalMass\t(3)FuelMass\t(4)OxidMass\t(5)CoM_x\t"
          "(6)J_x\t(7)J_y\t(8)JDot_x\t(9)JDot_y\t"
       << endl;
  cout << '#' << endl;

  //-------------------------------------------------------------------------//
  // Stage2 Params as a function of Flight Time:                             //
  //-------------------------------------------------------------------------//
  // XXX: Currently assuming no Thrust Vector Control, and Sea-Level Atmospheric
  // Pressure:
  S2::VernDeflections const vernDefls0; // All 0s by default

  // To verify the validity of Mass and MoI "Dots", we integrate them and compa-
  // re the results with the analytical Mass and MoIs:
  Mass   integrMass;
  MoI    integrMoIX;
  MoI    integrMoIY;
  Len    integrCoMX;
  double maxRelErrMass = 0.0;
  double maxRelErrMoIX = 0.0;
  double maxRelErrMoIY = 0.0;
  double maxRelErrCoMX = 0.0;

  constexpr Time   tau = 0.1_sec;

  for (Time t = 0.0_sec; t <= 300.0_sec; t += tau)
  {
    StageDynParams<LVSC::Soyuz21b> dp =
      S2::GetDynParams(SC::LiftOffTime + t, p0, vernDefls0);

    assert(IsZero(dp.m_com[1]) && IsZero(dp.m_com[2]) &&
           dp.m_mois      [1]  == dp.m_mois      [2]  &&
           dp.m_moiDots   [1]  == dp.m_moiDots   [2]);

    if (IsZero(t))
    {
      integrMass   = dp.m_fullMass;
      integrMoIX   = dp.m_mois[0];
      integrMoIY   = dp.m_mois[1];
      integrCoMX   = dp.m_com [0];
    }
    else
    {
      // Update the integrals using the backward-looking "Dots":
      integrMass  += tau * dp.m_fullMassDot;
      integrMoIX  += tau * dp.m_moiDots[0];
      integrMoIY  += tau * dp.m_moiDots[1];
      integrCoMX  += tau * dp.m_comDots[0];

      // Compute the integration errors:
      double integrMassErr = Abs(double(integrMass / dp.m_fullMass) - 1.0);
      double integrMoIErrX = Abs(double(integrMoIX / dp.m_mois[0] ) - 1.0);
      double integrMoIErrY = Abs(double(integrMoIY / dp.m_mois[1] ) - 1.0);
      double integrCoMErrX = Abs(double(integrCoMX / dp.m_com [0] ) - 1.0);

      maxRelErrMass = std::max(maxRelErrMass, integrMassErr);
      maxRelErrMoIX = std::max(maxRelErrMoIX, integrMoIErrX);
      maxRelErrMoIY = std::max(maxRelErrMoIY, integrMoIErrY);
      maxRelErrCoMX = std::max(maxRelErrCoMX, integrCoMErrX);
    }

    cout << t.Magnitude()               << '\t'
         << dp.m_fullMass  .Magnitude() << '\t'
         << dp.m_fuelMass  .Magnitude() << '\t'
         << dp.m_oxidMass  .Magnitude() << '\t'
         << dp.m_com    [0].Magnitude() << '\t'
         << dp.m_mois   [0].Magnitude() << '\t'
         << dp.m_mois   [1].Magnitude() << '\t'
         << dp.m_moiDots[0].Magnitude() << '\t'
         << dp.m_moiDots[1].Magnitude() << endl;
  }
  cout << "# Stage2FullMR        : " << S2::FullMR     << endl;
  cout << "# Stage2MinEndMass    : " << S2::MinEndMass << endl;
  cout << "# MaxRelErrMass       : " << maxRelErrMass  << endl;
  cout << "# MaxRelErrMoIX       : " << maxRelErrMoIX  << endl;
  cout << "# MaxRelErrMoIY       : " << maxRelErrMoIY  << endl;
  cout << "# MaxRelErrCoMX       : " << maxRelErrCoMX  << endl;
  return 0;
}
