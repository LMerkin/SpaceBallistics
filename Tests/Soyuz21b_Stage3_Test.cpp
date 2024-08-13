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
  cout << "# Stage3PropMargin  : "
       << (double(S3::MaxBurnDur / SC::Stage3BurnDur) - 1.0) * 100.0
       << " %" << endl;

  // Percentage of the StaticThrust in the over-all thrust of Stage3:
  cout << "# Stage3StaticThrust: " << S3::StaticThrust
       << "\n#\t\t\t(" << (double(S3::StaticThrust/S3::ThrustVac) * 100.0)
       << " % of over-all)" << endl;

  // Header for the following table:
  cout << "# Time\tTotalMass\tFuelMass\tOxidMass\tCoM_x\tJ_x\tJ_y" << endl;

  //-------------------------------------------------------------------------//
  // Stage3 Params as a function of Flight Time:                             //
  //-------------------------------------------------------------------------//
  // XXX: Currently assuming no Thrust Vector Control:
  //
  S3::ChamberDeflections chamberDefls0;    // All 0s by default 

  for (Time t = 250.0_sec; t <= 600.0_sec; t += 0.1_sec)
  {
    StageDynParams<LVSC::Soyuz21b> dp = S3::GetDynParams(t, chamberDefls0);

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
  return 0;
}
