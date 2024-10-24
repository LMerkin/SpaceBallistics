// vim:ts=2:et
//===========================================================================//
//                         "Tests/DE440T_Test.cpp":                          //
//===========================================================================//
#include <iostream>

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
  TDB  from  = DE440T::Bits::From;
  TDB  to    = DE440T::Bits::To;
  Time step  = To_Time(1.0_day);

  for (TDB tdb = from; tdb <= to; tdb += step)
  {
    TT   tt    = DE440T::TTofTDB(tdb);
    Time delta = tt.GetTime()  - tdb.GetTime();

    cout << (tdb - from).Magnitude() << '\t' << delta.Magnitude() << endl;
  }
  return 0;
}
