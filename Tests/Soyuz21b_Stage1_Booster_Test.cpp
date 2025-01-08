// vim:ts=2:et
//===========================================================================//
//                "Tests/Soyuz21b_Stage1_Booster_Test.cpp":                  //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"
#include "SpaceBallistics/Utils.hpp"
#include "TestUtils.hpp"

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
  cout << "# BoosterOxidTankLowD   : " << B::OxidTankLowD             << endl;
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
          "(6)CoM_y\t(7)Com_z\t(8)J_x\t(9)J_y\t(10)J_z\t(11)JDot_x\t"
          "(12)JDot_y\t(13)JDot_z" << endl;
  cout << '#' << endl;

  //-------------------------------------------------------------------------//
  // Booster Block Params as a function of Flight Time:                      //
  //-------------------------------------------------------------------------//
  // We need a well-defined LiftOff Time:
  FlightTime const LiftOffTime(0.0_sec, TTofStr("2025-01-01_00:00:00"));

  // XXX: Currently assuming no Thrust Vector Control, Sea-Level Atmospheric
  // Pressure, and no AeroDynamic Controls:
  //
  constexpr B::VernDeflections vernDefls0; // All 0s by default
  constexpr ME::VelVE          v0;         //
  constexpr ME::VelVE          wind0;      //
  constexpr Angle_deg          fin0;       //

  // To verify the validity of Mass and MoI "Dots", we integrate them and compa-
  // re the results with the analytical Mass and MoIs:
  Mass   integrMass;
  MoI    integrMoIX, integrMoIY, integrMoIZ;
  Len    integrCoMX, integrCoMY, integrCoMZ;

  double maxRelErrMass = 0.0;
  double maxRelErrMoIX = 0.0;
  double maxRelErrMoIY = 0.0;
  double maxRelErrMoIZ = 0.0;
  double maxRelErrCoMX = 0.0;
  double maxRelErrCoMY = 0.0;
  double maxRelErrCoMZ = 0.0;

  constexpr Time   tau = 0.1_sec;

  for (Time t = 0.0_sec; t <= 120.0_sec; t += tau)
  {
    StageDynParams<LVSC::Soyuz21b> dp =
      B::GetDynParams(LiftOffTime + t, p0, vernDefls0, v0, wind0, fin0);

    if (IsZero(t))
    {
      integrMass   = dp.m_fullMass;
      integrMoIX   = dp.m_mois[0];
      integrMoIY   = dp.m_mois[1];
      integrMoIZ   = dp.m_mois[2];
      integrCoMX   = dp.m_com [0];
      integrCoMY   = dp.m_com [1];
      integrCoMZ   = dp.m_com [2];
    }
    else
    {
      // Update the integrals using the backward-looking "Dots":
      integrMass  += tau * dp.m_fullMassDot;
      integrMoIX  += tau * dp.m_moiDots[0];
      integrMoIY  += tau * dp.m_moiDots[1];
      integrMoIZ  += tau * dp.m_moiDots[2];
      integrCoMX  += tau * dp.m_comDots[0];
      integrCoMY  += tau * dp.m_comDots[1];
      integrCoMZ  += tau * dp.m_comDots[2];

      // Compute the integration errors:
      double integrMassErr = Abs(double(integrMass / dp.m_fullMass) - 1.0);
      double integrMoIErrX = Abs(double(integrMoIX / dp.m_mois[0] ) - 1.0);
      double integrMoIErrY = Abs(double(integrMoIY / dp.m_mois[1] ) - 1.0);
      double integrMoIErrZ = Abs(double(integrMoIZ / dp.m_mois[2] ) - 1.0);
      double integrCoMErrX = Abs(double(integrCoMX / dp.m_com [0] ) - 1.0);
      double integrCoMErrY = Abs(double(integrCoMY / dp.m_com [1] ) - 1.0);
      double integrCoMErrZ = Abs(double(integrCoMZ / dp.m_com [2] ) - 1.0);

      maxRelErrMass = std::max(maxRelErrMass, integrMassErr);
      maxRelErrMoIX = std::max(maxRelErrMoIX, integrMoIErrX);
      maxRelErrMoIY = std::max(maxRelErrMoIY, integrMoIErrY);
      maxRelErrMoIZ = std::max(maxRelErrMoIZ, integrMoIErrZ);
      maxRelErrCoMX = std::max(maxRelErrCoMX, integrCoMErrX);
      maxRelErrCoMY = std::max(maxRelErrCoMY, integrCoMErrY);
      maxRelErrCoMZ = std::max(maxRelErrCoMZ, integrCoMErrZ);
    }

    cout << t.Magnitude()               << '\t'
         << dp.m_fullMass  .Magnitude() << '\t'
         << dp.m_fuelMass  .Magnitude() << '\t'
         << dp.m_oxidMass  .Magnitude() << '\t'
         << dp.m_com    [0].Magnitude() << '\t'
         << dp.m_com    [1].Magnitude() << '\t'
         << dp.m_com    [2].Magnitude() << '\t'
         << dp.m_mois   [0].Magnitude() << '\t'
         << dp.m_mois   [1].Magnitude() << '\t'
         << dp.m_mois   [2].Magnitude() << '\t'
         << dp.m_moiDots[0].Magnitude() << '\t'
         << dp.m_moiDots[1].Magnitude() << '\t'
         << dp.m_moiDots[2].Magnitude() << endl;
  }
  cout << "# BoosterFullMR         : " << B::FullMR      << endl;
  cout << "# BoosterMinEndMass     : " << B::MinEndMass  << endl;
  cout << "# MaxRelErrMass         : " << maxRelErrMass  << endl;
  cout << "# MaxRelErrMoIX         : " << maxRelErrMoIX  << endl;
  cout << "# MaxRelErrMoIY         : " << maxRelErrMoIY  << endl;
  cout << "# MaxRelErrMoIZ         : " << maxRelErrMoIZ  << endl;
  cout << "# MaxRelErrCoMX         : " << maxRelErrCoMX  << endl;
  cout << "# MaxRelErrCoMY         : " << maxRelErrCoMY  << endl;
  cout << "# MaxRelErrCoMZ         : " << maxRelErrCoMZ  << endl;
  return 0;
}
