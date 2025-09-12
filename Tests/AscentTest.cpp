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

  if (argc < 5)
  {
    cerr << "PARAMETERS: PayLoad_kg Perigee_km Apogee_km NParams"
         << endl;
    return 1;
  }
  Mass      mL      { atof         (argv[1])  };
  LenK      hq      { atof         (argv[2])  };
  LenK      hQ      { atof         (argv[3])  };
  unsigned  NP      { unsigned(atoi(argv[4])) };

  constexpr double    alpha1     = 4.11;
  constexpr ForceK    thrust2Vac = 1.0 * 63'700.0_kg * g0K;
  constexpr ForceK    thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  constexpr Angle_deg incl       = 64.0_deg;
  constexpr Angle_deg launchLat  = 63.0_deg;
  constexpr int       logLevel   = 0;
  constexpr int       maxEvals   = 10000;
  try
  {
    optional<Ascent2::AscCtlsD> optCtls =
      Ascent2::FindOptimalAscentCtls
      (
        alpha1, thrust2Vac, thrust1Vac,
        mL, hq, hQ,   incl, launchLat,
        &cout,  logLevel,
        Ascent2::OptMethod::NOMAD4, NP, maxEvals
      );

    if (!bool(optCtls))
      cout << "FIG!" << endl;
    else
    {
      Ascent2::AscCtlsD const& v = optCtls.value();
      cout
        << "T2        = " << v.m_T2        << endl
        << "aMu2      = " << v.m_aMu2      << endl
        << "bMu2      = " << v.m_bMu2      << endl
        << "upRate2   = " << v.m_upRate2   << endl
        << "aAoA2     = " << v.m_aAoA2     << endl
        << "bAoA2     = " << v.m_bAoA2     << endl << endl
        << "TGap      = " << v.m_TGap      << endl << endl
        << "T1        = " << v.m_T1        << endl
        << "aMu1      = " << v.m_aMu1      << endl
        << "bMu1      = " << v.m_bMu1      << endl
        << "upRate1   = " << v.m_upRate1   << endl
        << "aAoA1     = " << v.m_aAoA1     << endl
        << "bAoA1     = " << v.m_bAoA1     << endl << endl
        << "startH    = " << v.m_startH    << endl
        << "startV    = " << v.m_startV    << endl
        << "startMass = " << v.m_startMass << endl
        << "ObjVal    = " << v.m_objVal    << endl;
    }
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
