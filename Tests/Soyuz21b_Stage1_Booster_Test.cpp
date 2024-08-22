// vim:ts=2:et
//===========================================================================//
//                "Tests/Soyuz21b_Stage1_Booster_Test.cpp":                  //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  using       B = Soyuz21b_Stage1_Booster<'B'>;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  cout << "# Stage1MaxFullThrDur: " << B::MaxFullThrustDur << endl;
  cout << "# Stage1CutOffTime   : " << B::CutOffTime       << endl;
  cout << "# Stage1MaxBurnTime  : " << B::MaxBurnTime      << endl;

  // Header for the following table:
  cout << '#' << endl;
  cout << "# Time\tTotalMass\tFuelMass\tOxidMass\tCoM_x\tJ_x\tJ_y" << endl;
  cout << '#' << endl;
  return 0;
}
