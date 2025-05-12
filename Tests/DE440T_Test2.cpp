// vim:ts=2:et
//===========================================================================//
//                        "Tests/DE440T_Test2.cpp":                          //
//       Equatorial (J2000.0) Ephemerides of the Sun in (alpha, delta)       //
//===========================================================================//
// Asserts must be enabled for this test:
#ifdef   NDEBUG
#undef   NDEBUG
#endif
#include "SpaceBallistics/CoOrds/TimeScales.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.hpp"
#include "SpaceBallistics/CoOrds/GeoCDynEqFixCOS.h"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"
#include "SpaceBallistics/PhysEffects/DE440T.h"
#include "SpaceBallistics/PhysEffects/EarthRotationModel.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include "SpaceBallistics/AngleUtils.hpp"
#include "TestUtils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>

using namespace SpaceBallistics;
using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2 || argc > 5)
  {
    cerr << "PARAMS: [-d] TT_From [TT_To [TimeStep_days]]" << endl;
    return 1;
  }
  // Use Dynamic Equator and (Mean) Ecliptic instead of GCRS?
  bool useDyn = strcmp(argv[1], "-d") == 0;
  int  off    = int(useDyn);

  // Get the TT range and step:
  TT   from   =                     TTofStr  (argv[off + 1]);
  TT   to     = (argc >= off + 3) ? TTofStr  (argv[off + 2])
                                  : from;
  Time step   = (argc >= off + 4) ? TimeOfStr(argv[off + 3])
                                  : To_Time(Time_day(1.0));
  if (from == TT() || to == TT() || from > to || !IsPos(step))
  {
    cerr << "Invalid TT Range or Step" << endl;
    return 1;
  }

  // Create the "GeoCDynEqFixCOS"s (only used if "useDyn" is set); they are not
  // assignable, so we may need to create multiple copies of them:
  vector<GeoCDynEqFixCOS> dynCOSes;
  if (useDyn)
    dynCOSes.emplace_back(from.GetJYr());

  // Generate the Ephemerides of the Sun in the J2000.0 Equatorial Co-Ords:
  for (TT tt = from; tt <= to; tt += step)
  {
    Time_jyr jyr = tt.GetJYr();
    TDB tdb (tt);

    PosKV_BCRS<Body::Earth> posE;
    PosKV_BCRS<Body::Sun>   posS;
    VelKV_BCRS<Body::Earth> velE;
    VelKV_BCRS<Body::Sun>   velS;
    DE440T::GetPlanetBarEqPV<Body::Earth>(tdb, &posE, &velE);
    DE440T::GetPlanetBarEqPV<Body::Sun>  (tdb, &posS, &velS);

    // Compute the GeoEq PV vectors of the Sun:
    PosKV_GCRS<Body::Sun>   posES = posS - posE;
    VelKV_GCRS<Body::Sun>   velES = velS - velE;

    // For the results (initially all 0s): Spgerical CoOrds in HMS and DMS:
    std::tuple<        Angle_hh,  Angle_mm,     Angle_ss>     hms;
    std::tuple<double, Angle_deg, Angle_arcMin, Angle_arcSec> dms;
    AngVel    alphaDot, deltaDot;
    VelK      radVel;

    if (useDyn)
    {
      if (Abs(jyr - dynCOSes.back().GetEpoch()) > 1.0_jyr)
        // Create a "more modern" "dynCOS":
        dynCOSes.emplace_back(jyr);

      PosKV_GeoDynEqFix<Body::Sun> posDynES =
        dynCOSes.back().ToDynEquinox(posES);

      VelKV_GeoDynEqFix<Body::Sun> velDynES =
        dynCOSes.back().ToDynEquinox(velES);

      SpherPV<GeoCDynEqFixCOS, Body::Sun> spherDynPVS(posDynES, velDynES);

      hms      = ToHMS(spherDynPVS.GetAlpha());
      dms      = ToDMS(spherDynPVS.GetDelta());
      alphaDot = spherDynPVS.GetAlphaDot();
      deltaDot = spherDynPVS.GetDeltaDot();
      radVel   = spherDynPVS.GetRadVel  ();
    }
    else
    {
      // Compute the GCRS Spherical CoOrds:
      GeoCEqSpherPV<Body::Sun> spherPVS(posES, velES);

      hms      = ToHMS(spherPVS.GetAlpha());
      dms      = ToDMS(spherPVS.GetDelta());
      alphaDot = spherPVS.GetAlphaDot();
      deltaDot = spherPVS.GetDeltaDot();
      radVel   = spherPVS.GetRadVel  ();

      // CHECK: Get back to the GeoEq PV vectors:
      // XXX: Although we explicitly reset the NDEBUG flag at the beginning, in
      // CLang it may still be set (???), so use the following guard to prevent
      // compile-time warnings.   However, this will prevent the following test
      // from being done:
#     ifndef NDEBUG
      auto  [posES1, velES1] = spherPVS.GetPVVectors();
      assert(posES1.x().ApproxEquals(posES.x()) &&
             posES1.y().ApproxEquals(posES.y()) &&
             posES1.z().ApproxEquals(posES.z()) &&
             velES1.x().ApproxEquals(velES.x()) &&
             velES1.y().ApproxEquals(velES.y()) &&
             velES1.z().ApproxEquals(velES.z()));
#     endif
    }
    // Compute the RA and Decl "Dots" (in Seconds / ArcSeconds per Day) in
    // GCRS or Dyn CoOrds:
    auto alphaDotSD = To_Time_day(To_Angle_ss    (alphaDot));
    auto deltaDotSD = To_Time_day(To_Angle_arcSec(deltaDot));

    // Get the UTC from this TT:
    UTC utc(tt);

    // Output:
    printf("%04d-%02d-%02d_%02d:%02d:%03.3lf\t"
             "%02.0lf:%02.0lf:%02.06lf\t"
           "%c%02.0lf_%02.0lf_%02.06lf\t"
             "%.3lf\t%.2lf\t%.2lf\t%.3lf\n",
           utc.m_year, utc.m_month, utc.m_day,
           utc.m_hour, utc.m_min,   utc.m_sec,
           get<0>(hms).Magnitude(), get<1>(hms).Magnitude(),
           get<2>(hms).Magnitude(),
           get<0>(dms) < 0  ? '-' : get<0>(dms) > 0 ? '+' : ' ',
           get<1>(dms).Magnitude(), get<2>(dms).Magnitude(),
           get<3>(dms).Magnitude(),
           LenK(posES).Magnitude(),
           alphaDotSD.Magnitude(),  deltaDotSD.Magnitude(),
           radVel.Magnitude());
  }
  // All Done!
  return 0;
}
