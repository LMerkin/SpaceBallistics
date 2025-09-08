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

  if (argc < 3)
  {
    cerr << "PARAMETERS: LEO_Altitude_km PayLoad_kg" << endl;
    return 1;
  }
  LenK   hC { atof(argv[1]) };
  Mass   mL { atof(argv[2]) };

  constexpr ForceK thrust2Vac = 1.0 * 63'700.0_kg * g0K;
  constexpr ForceK thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  constexpr double alpha1     = 4.11;
  try
  {
    Ascent2 asc
    (
      alpha1, thrust2Vac,   thrust1Vac,
      mL, hC, hC, 64.0_deg, 63.0_deg, Ascent2::AscCtls{},
      &cout,  0
    );
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
