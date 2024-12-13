// vim:ts=2:et
//===========================================================================//
//                         "Tests/TT_UTC_Test.cpp":                          //
//===========================================================================//
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <cstdlib>
#include <cstdio>

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  if (argc != 1 && argc != 2)
  {
    fputs("PARAM: [secFrac in [0..1)]\n", stderr);
    return 1;
  }

  //-------------------------------------------------------------------------//
  // UTC->TAI Conversion Test:                                               //
  //-------------------------------------------------------------------------//
  for (int year  = 1960; year  <= 2025; ++year)
  for (int month = 1;    month <= 12;   ++month)
  {
    UTC utc(year, month, 1,  0, 0, 0.0);
    auto [jd_utc, DeltaAT] = TBits::MkTAI(utc);

    double yf  = double(year) + double(month - 1)/12.0;

    printf("%.3lf\t%.2lf\t%.3lf\n",
           yf, jd_utc.Magnitude(), DeltaAT.Magnitude());

    // Convert this UTC to TT, back to UTC and again to TT (since we don't want
    // to compare the UTCs due to day/month/year wrappings).
    // XXX: There is a "DeltaAT" gap @ 1961-01-01, so do not perform the UTCs
    // comparison at that point, since the result may be affected by rounding
    // errors:
    if (!(year == 1961 && month == 1))
    {
      TT  tt0 (utc);
      UTC utc1(tt0);
      TT  tt1 (utc1);
      assert(tt0.ApproxEquals(tt1, 1e-10));
    }
  }

  //-------------------------------------------------------------------------//
  // Leap Seconds Test:                                                      //
  //-------------------------------------------------------------------------//
  double secFrac = (argc == 2) ? atof(argv[1]) : 0.0;
  if (secFrac < 0 || secFrac >= 1)
  {
    printf("ERROR: Invalid secFrac=%lf\n", secFrac);
    return 1;
  }

  for (int i = 0; i < UTC::NLS; ++i)
  {
    auto [y, m, d] = UTC::LeapSecondDates[i];

    // It may only be end of June or end of December:
    assert((m == 6 && d == 30) || (m == 12 && d == 31));

    // Into this Leap Second:
    UTC utc0(y, m, d, 23, 59, 60.0 + secFrac);

    // Into the next second:
    UTC utc1
    (
      (m == 6) ? y : (y+1),
      (m == 6) ? 7 : 1,
      1,
      0, 0, secFrac
    );
    // Compute the TT for both instants:
    TT tt0(utc0);
    TT tt1(utc1);

    // They should be separated by (almost) exactly 1 second:
    DEBUG_ONLY(Time   diff = tt1 - tt0;)
    assert(diff.ApproxEquals(1.0_sec));

    // Convert the both to UTC and back:
    UTC utc2(tt0);
    UTC utc3(tt1);
    TT  tt2 (utc2);
    TT  tt3 (utc3);
    assert( tt2.ApproxEquals(tt0));
    assert( tt3.ApproxEquals(tt1));
    assert((tt3 - tt2).ApproxEquals(1.0_sec));
  }
}
