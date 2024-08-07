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

  cout << "# Stage2PropMargin  : "
       << (double(S2::MaxFlightTime / SC::Stage2CutOffTime) - 1.0) * 100.0
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

  // Over-All Length:
  cout << "# Stage2OverAllLen  : " << S2::H << endl;

  cout << "# StillFree         : "
       << (S2::FuelTankBtm.GetLow()[0] - (S2::TopX - S2::H)) << endl;

  cout << "# FromTopToTail     : "
       << (S2::TailEncl.GetLow()[0]    - (S2::TopX - S2::H)) << endl;
  return 0;
}
