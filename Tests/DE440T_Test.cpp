// vim:ts=2:et
//===========================================================================//
//                         "Tests/DE440T_Test.cpp":                          //
//===========================================================================//
#include <cstdio>

// Asserts must be enabled for this test:
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "SpaceBallistics/PhysForces/DE440T.hpp"

int main()
{
  using namespace SpaceBallistics;
  using namespace std;

  // Self-Consistency Test:
  DE440T::SelfTest();

  // TT-TDB plot:
  Time step  = To_Time(0.1_day);
  TDB  from  = DE440T::Bits::From;
  TDB  to    = DE440T::Bits::To - step;  // Safety margin

  for (TDB tdb = from; tdb <= to; tdb += step)
  {
    TT   tt    = DE440T::TTofTDB(tdb);
    Time delta = tt.GetTime()  - tdb.GetTime();

    // Compute and verify the inverse.
    // We should achieve at least a nanosecond precision:
    TDB  tdb1  = DE440T::TDBofTT(tt);
    Time err   = Abs(tdb1 - tdb);
    assert(err < Time(1e-9));

    printf("%.1lf\t%.6lf\t%.6lf\t%.1lf\t%.6lf\t%.6lf\n",
           (tdb - from).Magnitude(),  delta.Magnitude(),
           err.Magnitude(),          tdb .GetTime().Magnitude(),
           tt.GetTime().Magnitude(), tdb1.GetTime().Magnitude());
  }
  return 0;
}
