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
/*
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
  constexpr int       maxEvals   = 10000;
  try
  {
    optional<pair<Ascent2::AscCtlsL, Ascent2::AscCtlsD>> optCtls =
      Ascent2::FindOptimalAscentCtls
      (
        alpha1, thrust2Vac, thrust1Vac,
        mL, hq, hQ,   incl, launchLat,
        &cout,  0,
        NP,     maxEvals
      );

    if (!bool(optCtls))
    {
      cout << "FIG!" << endl;
      return 1;
    }

    // If OK:
    // Output the optimal params:
    //
    Ascent2::AscCtlsD const& v = optCtls.value().second;
    cout
      << "T2        = " << v.m_T2         << endl
      << "aMu2      = " << v.m_aMu2       << endl
      << "bMu2      = " << v.m_bMu2       << endl
      << "upRate2   = " << v.m_upRate2    << endl
      << "aAoA2     = " << v.m_aAoA2      << endl
      << "bAoA2     = " << v.m_bAoA2      << endl << endl
      << "TGap      = " << v.m_TGap       << endl << endl
      << "T1        = " << v.m_T1         << endl
      << "aMu1      = " << v.m_aMu1       << endl
      << "bMu1      = " << v.m_bMu1       << endl
      << "upRate1   = " << v.m_upRate1    << endl
      << "aAoA1     = " << v.m_aAoA1      << endl
      << "bAoA1     = " << v.m_bAoA1      << endl << endl;

    // Now make the final run with optimal params found and produce the output:
    Ascent2 asc
    (
      alpha1, thrust2Vac, thrust1Vac,
      mL, hq, hQ,   incl, launchLat,
      &cout,  3
    );
    asc.SetCtlParams(optCtls.value().first);
    Ascent2::RunRes res = asc.Run();

    cout
      << "RESULT    = " << Ascent2::ToString(res.m_rc) << endl
      << "startTime = " << res.m_T        << endl
      << "startH    = " << res.m_hT       << endl
      << "startV    = " << res.m_VT       << endl
      << "startMass = " << res.m_mT       << endl
      << "maxQ      = " << res.m_maxQ     << endl
      << "sepQ      = " << res.m_sepQ     << endl
      << "maxLongG  = " << res.m_maxLongG << endl << endl;
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
*/
	return  0;
}
