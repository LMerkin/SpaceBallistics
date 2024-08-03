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

  // Max and Real Fuel and Oxid Loads:
  cout << "# Stage2MaxOxidLoad : " << S2::OxidTankMC             << endl;
  cout << "# Stage2ActOxidLoad : " << S2::OxidLoadRatio * 100.0  << " %"
       << endl;

  return 0;
}
