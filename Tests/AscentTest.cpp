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
    cerr << "PARAMETERS: PayLoad_kg Perigee_km Apogee_km" << endl;
    return 1;
  }
  Mass   mL { atof(argv[1]) };
  LenK   hq { atof(argv[2]) };
  LenK   hQ { atof(argv[3]) };

  constexpr double    alpha1     = 4.11;
  constexpr ForceK    thrust2Vac = 1.0 * 63'700.0_kg * g0K;
  constexpr ForceK    thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  constexpr Angle_deg incl       = 64.0_deg;
  constexpr Angle_deg launchLat  = 63.0_deg;
  constexpr int       logLevel   = 0;
  try
  {
    std::optional<Ascent2::AscCtlsD> optCtls =
      Ascent2::FindOptimalAscentCtls
      (
        alpha1, thrust2Vac, thrust1Vac,
        mL, hq, hQ,   incl, launchLat,
        &cout,  logLevel
      );

    if (!bool(optCtls))
      std::cout << "FIG!" << std::endl;
    else
    {
      Ascent2::AscCtlsD const& v = optCtls.value();
      std::cout
        << v.m_T2       .Magnitude() << "  "
        << v.m_aMu2     .Magnitude() << "  "
        << v.m_bMu2     .Magnitude() << "  "
        << v.m_aAoA2    .Magnitude() << "  "
        << v.m_bAoA2    .Magnitude() << "  "
        << v.m_TGap     .Magnitude() << "  "
        << v.m_T1       .Magnitude() << "  "
        << v.m_aMu1     .Magnitude() << "  "
        << v.m_bMu1     .Magnitude() << "  "
        << v.m_aAoA1    .Magnitude() << "  "
        << v.m_bAoA1    .Magnitude() << "  "
        << v.m_startMass.Magnitude() << std::endl;
    }
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
