// vim:ts=2:et
//===========================================================================//
//                     "Tests/Soyuz21b_Booster_Test.cpp":                    //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Booster.h"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  using       B = Soyuz21b_Booster<'B'>;

  //-------------------------------------------------------------------------//
  // Geometry and Mass Params:                                               //
  //-------------------------------------------------------------------------//
  cout << "# Stage1MaxFullThrTime: " << B::MaxFullThrustTime      << endl;
  cout << "# Stage1CutOffTime    : " << B::CutOffTime             << endl;
  cout << "# Stage1MaxBurnTime   : " << B::MaxBurnTime            << endl;

  // Header for the following table:
  cout << "# Time\tTotalMass\tFuelMass\tOxidMass\tCoM_x\tJ_x\tJ_y" << endl;
  return 0;
}
