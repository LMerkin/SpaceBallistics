// vim:ts=2:et
//===========================================================================//
//                             "Stage2Test.cpp":                             //
//===========================================================================//
#include "SpaceBallistics/LVSC/Soyuz21b_Stage2.h"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;
  using namespace SpaceBallistics;
  namespace  SC = Soyuz21b_Consts;
  using      S2 = Soyuz21b_Stage2;

  cout << "Stage2PropMargin: "
       << (double(S2::MaxFlightTime / SC::Stage2CutOffTime) - 1.0) * 100.0
       << " %" << endl;
  return 0;
}
