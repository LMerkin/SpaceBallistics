// vim:ts=2:et
//===========================================================================//
//                         "Tests/DE440T_Test0.cpp":                         //
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

  Body body   = Body::Sun;    // Just to initialise...
  bool isNut  = false;
  bool isLibr = false;
  bool isTT   = false;

  if (strcmp(argv[1], "EarthNutations") == 0)
    isNut  = true;
  else
  if (strcmp(argv[1], "MoonLibrations") == 0)
    isLibr = true;
  else
  if (strcmp(argv[1], "TT_TDB") == 0)
    isTT = true;
  else
  try
  {
    // Get the actual Body:
    body = ToBody(argv[1]);
  }
  catch (...)
  {
    cerr << "ERROR: Invalid Body: " << argv[1] << endl;
    return 1;
  }

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
    if (isNut)
    {
      Angle   nuts[2];
      DE440T::GetEarthNutations(tdb, nuts);

      printf("%.2lf\t%.16e\t%.16e\n",
             tdb.GetJD().Magnitude(), double(nuts[0]), double(nuts[1]));
    }
    else
    if (isLibr)
    {
      Angle   librs   [3];
      AngVel  librDots[3];
      DE440T::GetMoonLibrations(tdb, librs, librDots);

      printf("%.2lf\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             tdb.GetJD().Magnitude(),
             double(librs[0]),        double(librs[1]),
             double(librs[2]),
             librDots[0].Magnitude(), librDots[1].Magnitude(),
             librDots[2].Magnitude());
    }
    else
    if (isTT)
    {
      TT   tt(tdb);
      Time delta = tt.GetTime()  - tdb.GetTime();

      // Compute and verify the inverse.
      // We should achieve at least a nanosecond precision:
      TDB  tdb1(tt);
      Time err   = Abs(tdb1 - tdb);
      assert(err < Time(1e-9));

      printf("%.2lf\t%.1lf\t%.6lf\t%.6lf\t%.1lf\t%.6lf\t%.6lf\n",
             tdb.GetJD() .Magnitude(),
             (tdb - from).Magnitude(),  delta.Magnitude(),
             err.Magnitude(),          tdb .GetTime().Magnitude(),
             tt.GetTime().Magnitude(), tdb1.GetTime().Magnitude());
    }
    else
    if (body == Body::Moon)
    {
      // For the Moon, get the GeoCentric Fixed-Axes Ecliptical CoOrds:
      PosKVGeoEclFix pos;
      VelKVGeoEclFix vel;
      DE440T::GetMoonGEclPV(tdb, &pos, &vel);

      // Output:
      printf("%.2lf\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
             tdb.GetJD().Magnitude(),
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
      printf("%.2lf\t"
             "%.16e\t%.16e\t%.16e\t%.16e\t"
             "%.16e\t%.16e\t%.16e\t%.16e\n",
             tdb.GetJD().Magnitude(),
             pos.x().Magnitude(), pos.y().Magnitude(), pos.z().Magnitude(),
             LenK(pos).Magnitude(),
             vel.x().Magnitude(), vel.y().Magnitude(), vel.z().Magnitude(),
             VelK(vel).Magnitude());
    }
  // All Done!
  return 0;
}
