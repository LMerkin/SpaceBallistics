// vim:ts=2:et
//===========================================================================//
//                         "Tests/DE440T_Test.cpp":                          //
//===========================================================================//
#include <cstdio>
#include <cstdlib>

// Asserts must be enabled for this test:
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "SpaceBallistics/PhysForces/DE440T.h"

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  //-------------------------------------------------------------------------//
  // Temporal Consistency Test:                                              //
  //-------------------------------------------------------------------------//
  DE440T::Bits::TemporalConsistencyTest();

  //-------------------------------------------------------------------------//
  // Test for a Particular Body:                                             //
  //-------------------------------------------------------------------------//
  if (argc != 3 && argc != 4 && argc != 5)
  {
    cerr << "PARAMS: BodyName TDB_From [TDB_To [TimeStep_days]]" << endl;
    return 1;
  }

  // Initialisation:
  Body body = Body::Sun;
  bool isTT = false;

  if (strcmp(argv[1], "TT_TDB") == 0)
    isTT = true;
  else
    // Get the actual Body:
    body = ToBody(argv[1]);

  // Get the JD_TDB:
  TDB  from    {Time_day(atof(argv[2]))};
  TDB  to    = (argc >= 4) ?     TDB{Time_day(atof(argv[3]))} : from;
  Time step  = To_Time((argc >= 5) ? Time_day(atof(argv[4]))  : 1.0_day);

  if (from < DE440T::Bits::From || to > DE440T::Bits::To || from > to ||
      !IsPos(step))
  {
    cerr << "Invalid TDB Range or Step" << endl;
    return 1;
  }

  // Generate the Ephemerides:
  for (TDB tdb = from; tdb <= to; tdb += step)
    if (isTT)
    {
      TT   tt(tdb);
      Time delta = tt.GetTime()  - tdb.GetTime();

      // Compute and verify the inverse.
      // We should achieve at least a nanosecond precision:
      TDB  tdb1(tt);
      Time err   = Abs(tdb1 - tdb);
      assert(err < Time(1e-9));

      printf("%.1lf\t%.6lf\t%.6lf\t%.1lf\t%.6lf\t%.6lf\n",
             (tdb - from).Magnitude(),  delta.Magnitude(),
             err.Magnitude(),          tdb .GetTime().Magnitude(),
             tt.GetTime().Magnitude(), tdb1.GetTime().Magnitude());
    }
    else
    if (body == Body::Moon)
    {
      // For the Moon, get the GeoCentrix Fixed-Axes Ecliptical CoOrds:
      PosKVGeoEclFix pos;
      VelKVGeoEclFix vel;
      DE440T::GetMoonGEclPV(tdb, &pos, &vel);

      // Output:
      printf("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             pos.x().Magnitude(), pos.y().Magnitude(), pos.z().Magnitude(),
             LenK(pos).Magnitude(),
             vel.x().Magnitude(), vel.y().Magnitude(), vel.z().Magnitude(),
             VelK(vel).Magnitude());
    }
    else
    {
      // Any other Body: It's Sun or a Planet, or EMB:
      // Get the BaryCentric Ecliptical CoOrds:
      PosKVBEcl pos;
      VelKVBEcl vel;
      DE440T::GetPlanetBEclPV(body, tdb, &pos, &vel);

      // Output:
      printf("%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             pos.x().Magnitude(), pos.y().Magnitude(), pos.z().Magnitude(),
             LenK(pos).Magnitude(),
             vel.x().Magnitude(), vel.y().Magnitude(), vel.z().Magnitude(),
             VelK(vel).Magnitude());
    }
  // All Done!
  return 0;
}
