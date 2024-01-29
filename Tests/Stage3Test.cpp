// vim:ts=2:et
//===========================================================================//
//                            "Tests/Test0.cpp":                             //
//===========================================================================//
#include "SpaceBallistics/Soyuz21b_Stage3.h"
#include "SpaceBallistics/Soyuz21b_Stage3_Consts.h"
#include <iostream>

int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  namespace S3C = Soyuz21b_Stage3_Consts;

  //-------------------------------------------------------------------------//
  // The Stage3 Object, local to "main", created on stack:                   //
  //-------------------------------------------------------------------------//
  Soyuz21b_Stage3 st3;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  SpherSegm const& oUp  = st3.GetOxidTankUp ();
  SpherSegm const& fLow = st3.GetFuelTankLow();

  // Gap in Stage3 between the bottom of FuelTank and the top of OxidTank:
  Len gap = oUp.GetLeft()[0] - fLow.GetRight()[0];
  cout << "# Stage3FuelOxidGap   : " << gap << endl;

  // The max possible fuel and oxidiser load for Stage3 (using the max geometr-
  // ical tank volumes):
  cout << "# Stage3MaxFuelLoad   : " << st3.m_maxFuelMass << endl;
  cout << "# Stage3ActFuelLoad   : "
       << (double(S3C::FuelMass / st3.m_maxFuelMass) * 100.0) << " %"  << endl;

  cout << "# Stage3MaxOxidLoad   : " << st3.m_maxOxidMass << endl;
  cout << "# Stage3ActOxidLoad   : "
       << (double(S3C::OxidMass / st3.m_maxOxidMass) * 100.0) << " %"  << endl;

  cout << "# Stage3MassRateFT    : " << S3C::MassRateFT     << endl;
  cout << "# Stage3MaxBurnTimeFT : " << st3.m_maxBurnTimeFT << endl;

  // Percentage of the StaticThrust in the over-all thrust of Stage3:
  cout << "# Stage3StaticThrust  : " << S3C::StaticThrust
       << "\n\t\t#\t(" << (double(S3C::StaticThrust/S3C::ThrustVac) * 100.0)
       << " % of over-all)" << endl;

  //-------------------------------------------------------------------------//
  // Stage3 Mass as a function of Flight Time:                               //
  //-------------------------------------------------------------------------//
  for (Time t = 0.0_sec; t <= 600.0_sec; t += 0.1_sec)
  {
    StageDynParams dp = st3.GetDynParams(t);

    cout << t.Magnitude  ()  << '\t'  << dp.m_fullMass.Magnitude() << '\t'
         << dp.m_fuelMass.Magnitude() << '\t'
         << dp.m_oxidMass.Magnitude() << std::endl;

    // The following invariant must hold:
    if (t < st3.AftJetTime)
      assert(Abs(dp.m_fullMass  - dp.m_fuelMass - dp.m_oxidMass -
                 S3C::GasesMass - S3C::EmptyMass) < 1.0_kg);
    else
      assert(Abs(dp.m_fullMass  - dp.m_fuelMass - dp.m_oxidMass -
                 S3C::GasesMass - (S3C::EmptyMass - S3C::AftMass)) < 1.0_kg);
  }
  return 0;
}
