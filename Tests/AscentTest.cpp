// vim:ts=2:et
//===========================================================================//
//                          "Tests/AscentTest.cpp":                          //
//         Ascent to Low Earth Orbit Integration in Various Modes            //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  if (argc < 4)
  {
    cerr << "PARAMETERS: LEO_Altitude_km PayLoad_kg m1/m2_Ratio" << endl;
    return 1;
  }
  LenK   hC     { atof(argv[1]) };
  Mass   mL     { atof(argv[2]) };
  double alpha1 { atof(argv[3]) };

//ForceK thrust2Vac = 1.0 * 70'000.0_kg * g0K;
  ForceK thrust2Vac = 1.0 * 90'000.0_kg * g0K;
  ForceK thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  try
  {
    Ascent2 asc(mL, alpha1, thrust2Vac, thrust1Vac, &cout, 0);

    asc.FindOptimalAscentCtls(hC, hC, 64.0_deg, 63.0_deg);
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
