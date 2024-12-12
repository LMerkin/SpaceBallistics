// vim:ts=2:et
//===========================================================================//
//                           "Tests/JD_Test.cpp":                            //
//===========================================================================//
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <cstdio>

int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  //-------------------------------------------------------------------------//
  // UTC->TAI Conversion Test:                                               //
  //-------------------------------------------------------------------------//
  for (int year  = 1960; year  <= 2025; ++year)
  for (int month = 1;    month <= 12;   ++month)
  {
    UTC utc(year, month, 1,  0, 0, 0.0);
    auto [JD_UTC, DeltaAT] = TBits::MkTAI(utc);

    double yf  = double(year) + double(month - 1)/12.0;

    printf("%.3lf\t%.2lf\t%.3lf\n",
           yf, JD_UTC.Magnitude(), DeltaAT.Magnitude());
  }

  //-------------------------------------------------------------------------//
  // Leap Seconds Test:                                                      //
  //-------------------------------------------------------------------------//
  for (int i = 0; i < UTC::NLS; ++i)
  {
    auto [y, m, d] = UTC::LeapSecondDates[i];

    // It may only be end of June or end of December:
    assert((m == 6 && d == 30) || (m == 12 && d == 31));

    // Beginning of this Leap Second:
    UTC utc0(y, m, d, 23, 59, 60.0);

    // Beginning of the next second:
    UTC utc1
    (
      (m == 6) ? y : (y+1),
      (m == 6) ? 7 : 1,
      1,
      0, 0, 0.0
    );

    // Compute the TT for both instants:
    TT tt0(utc0);
    TT tt1(utc1);

    // They should be separated by (almost) exactly 1 second:
    DEBUG_ONLY(Time   diff = tt1 - tt0;)
    assert(diff.ApproxEquals(1.0_sec));
  }
}
