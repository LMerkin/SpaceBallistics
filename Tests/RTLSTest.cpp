// vim:ts=2:et
//===========================================================================//
//                          "Tests/RTLSTest.cpp":                            //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include <iostream>

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  if (argc != 2)
  {
    cerr << "PARAMETER: Config.ini" << endl;
    return 1;
  }
  string configIni(argv[1]);

  try
  {
    auto [optRes, finalRunRes] =
      RTLS1::FindOptimalReturnCtls(configIni, &cout);

    if (!bool(optRes))
    {
      // Optimisation Failed:
      cout << "FIG!" << endl;
      return 1;
    }

    // If OK: Output the optimal params:
    cout << "OPTIMAL SOLUTION FOUND:"       << endl;
    cout << optRes.value()                  << endl;

    if (bool(finalRunRes))
    {
      auto const& frr = finalRunRes.value();
      cout
        << "FINAL RUN:"                     << endl
        << "RESULT    = " << RTLS1::Base::ToString(frr.m_rc) << endl
        << "finalTime = " << frr.m_T        << endl
        << "finalL    = " << frr.m_LT       << endl
        << "finalV    = " << frr.m_VT       << endl
        << "finalMass = " << frr.m_mT       << endl
        << "maxQ      = " << frr.m_maxQ     << endl
        << "sepQ      = " << frr.m_sepQ     << endl
        << "maxLongG  = " << frr.m_maxLongG << endl << endl;
    }
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
