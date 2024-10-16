// vim:ts=2:et
//===========================================================================//
//                    "Tests/Soyuz21b_State3_Test.cpp":                      //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"
#include <iostream>

int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  namespace  SC = Soyuz21b_Consts;
  using      S3 = Soyuz21b_Stage3;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  cout << "# Stage3FuelOxidGap : " << S3::FuelOxidTanksGap << endl;

  // Max and Real Fuel and Oxid Loads:
  cout << "# Stage3MaxFuelLoad : " << S3::FuelTankMC             << endl;
  cout << "# Stage3ActFuelLoad : " << S3::FuelLoadRatio * 100.0  << " %"
       << endl;
  cout << "# Stage3MaxOxidLoad : " << S3::OxidTankMC             << endl;
  cout << "# Stage3ActOxidLoad : " << S3::OxidLoadRatio * 100.0  << " %"
       << endl;

  // Fuel and Oxid Masses at CutOff (as compared to UnSpendable Remnants):
  cout << "# Stage3CutOffFuel  : " << S3::CutOffFuelMass
       << "\n#\t\t\t(UnSpendableRemnant: " << S3::FuelRem << ')' << endl;
  cout << "# Stage3CutOffOxid  : " << S3::CutOffOxidMass
       << "\n#\t\t\t(UnSpendableRemnant: " << S3::OxidRem << ')' << endl;

  // Max Theoretical and Real Burn Durations:
  cout << "# Stage3MaxBurnDur  : " << S3::MaxBurnDur    << endl;
  cout << "# Stage3BurnDur     : " << SC::Stage3BurnDur << endl;
  cout << "# Stage3OxidFuelRat : " << S3::OxidFuelRat   << endl;
  cout << "# Stage3FuelMassRate: " << S3::FuelMR        << endl;
  cout << "# Stage3OxidMassRate: " << S3::OxidMR        << endl;

  // Geometry:
  cout << "# Stage3NozzlesExt  : " << S3::NozzlesExtH   << endl;
  cout << "# Stage3NozzlesLowX : " << S3::NozzlesLowX   << endl;

  // Header for the following table:
  cout << '#' << endl;
  cout << "# (1)Time\t(2)TotalMass\t(3)FuelMass\t(4)OxidMass\t(5)CoM_x\t"
          "(6)J_x\t(7)J_y\t(8)J_z\t(9)JDot_x\t(10)JDot_y\t(11)JDot_z"
       << endl;
  cout << '#' << endl;

  //-------------------------------------------------------------------------//
  // Stage3 Params as a function of Flight Time:                             //
  //-------------------------------------------------------------------------//
  // XXX: Currently assuming no Thrust Vector Control:
  //
  S3::ChamberDeflections const chamberDefls0; // All 0s by default

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

  constexpr Time   t0  = 250.0_sec;
  constexpr Time   t1  = 600.0_sec;
  constexpr Time   tau = 0.1_sec;
  bool      willJump   = true;

  for (Time t = t0; t <= t1; t += tau)
  {
    StageDynParams<LVSC::Soyuz21b> dp = S3::GetDynParams(t, chamberDefls0);

    assert(IsZero(dp.m_com[1]) && IsZero(dp.m_com[2]) &&
           dp.m_mois      [1]  == dp.m_mois      [2]  &&
           dp.m_moiDots   [1]  == dp.m_moiDots   [2]  &&
           IsZero(dp.m_comDots[1]) && IsZero(dp.m_comDots[2]));

    // XXX: Because the Mass, MoIs and CoM experience a jump, we have to re-init
    // our integration at that point:
    if (t == t0 || (willJump && t >= S3::AftJetTime))
    {
      integrMass   = dp.m_fullMass;
      integrMoIX   = dp.m_mois[0];
      integrMoIY   = dp.m_mois[1];
      integrCoMX   = dp.m_com [0];
      if (t >= S3::AftJetTime)
        willJump   = false;
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

    cout << t.Magnitude  ()             << '\t'
         << dp.m_fullMass  .Magnitude() << '\t'
         << dp.m_fuelMass  .Magnitude() << '\t'
         << dp.m_oxidMass  .Magnitude() << '\t'
         << dp.m_com    [0].Magnitude() << '\t'
         << dp.m_mois   [0].Magnitude() << '\t'
         << dp.m_mois   [1].Magnitude() << '\t'
         << dp.m_moiDots[0].Magnitude() << '\t'
         << dp.m_moiDots[1].Magnitude() << endl;
  }
  cout << "# MaxRelErrMass     : " << maxRelErrMass  << endl;
  cout << "# MaxRelErrMoIX     : " << maxRelErrMoIX  << endl;
  cout << "# MaxRelErrMoIY     : " << maxRelErrMoIY  << endl;
  cout << "# MaxRelErrCoMX     : " << maxRelErrCoMX  << endl;

  return 0;
}
