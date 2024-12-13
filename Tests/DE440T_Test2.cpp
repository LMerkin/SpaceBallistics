// vim:ts=2:et
//===========================================================================//
//                        "Tests/DE440T_Test2.cpp":                          //
//       Equatorial (J2000.0) Ephemerides of the Sun in (alpha, delta)       //
//===========================================================================//
// Asserts must be enabled for this test:
#ifdef   NDEBUG
#undef   NDEBUG
#endif
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include "SpaceBallistics/Utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <ctime>

using namespace SpaceBallistics;
using namespace std;

namespace
{
  //=========================================================================//
  // "TTofStr":                                                              //
  //=========================================================================//
  TT TTofStr(char const* a_str)
  {
    assert(a_str != nullptr && *a_str != '\0');

    // First of all, if "a_str" contains the "JD" prefix, we assume it is TT
    // in the JD format:
    if (a_str[0] == 'J' && a_str[1] == 'D')
      return TT(Time_day(atof(a_str + 2)));

    // Otherwise, it must be a UTC time-stamp in the format
    // "YYYY-MM-DD_hh:mm:ss":
    tm tmUTC;
    if (*strptime(a_str, "%Y-%m-%d_%H:%M:%S", &tmUTC) != '\0')
      return TT();

    UTC utc
    (
      tmUTC.tm_year + 1900,
      tmUTC.tm_mon  + 1,
      tmUTC.tm_mday,
      tmUTC.tm_hour,
      tmUTC.tm_min,
      double(tmUTC.tm_sec)
    );
    return TT(utc);
  }

  //=========================================================================//
  // "TimeOfStr":                                                            //
  //=========================================================================//
  Time TimeOfStr(char const* a_str)
  {
    assert (a_str != nullptr && *a_str != '\0');
    switch (a_str[0])
    {
      case 's': return Time            (atof(a_str + 1));
      case 'd': return To_Time(Time_day(atof(a_str + 1)));
      case 'y': return To_Time(Time_jyr(atof(a_str + 1)));
      default : return Time();
    }
  }
}

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

  // Generate the Ephemerides of the Sun in the J2000.0 Equatorial Co-Ords:
  for (TT tt = from; tt <= to; tt += step)
  {
    TDB tdb(tt);

    PosKVBarEq<Body::Earth> posE;
    PosKVBarEq<Body::Sun>   posS;
    VelKVBarEq<Body::Earth> velE;
    VelKVBarEq<Body::Sun>   velS;
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
      // Use "Dynamic Ecliptic and Mean Equator of the Date", so perform a man-
      // ual conversion of GCRS into DynEq (XXX: we don't use "GeoCDynEqFixCOS"
      // here yet): Extract the rotation matrices from the corresp ERM:
      EarthRotationModel   erm(tt);
      Mtx33 const& toDyn = erm.GetInvPN();

      // XXX: have to use COS=void in the corresp Vectors yet:
      Vector3D<LenK, void, Body::Sun> posDynES;
      Vector3D<VelK, void, Body::Sun> velDynES;

      MVMult33(toDyn, posES.GetArr(), posDynES.GetArr());
      MVMult33(toDyn, velES.GetArr(), velDynES.GetArr());

      SpherPV<void, Body::Sun> spherDynPVS(posDynES, velDynES);
      hms      = ToHMS(spherDynPVS.GetAlpha());
      dms      = ToDMS(spherDynPVS.GetDelta());
      alphaDot = spherDynPVS.GetAlphaDot();
      deltaDot = spherDynPVS.GetDeltaDot();
      radVel   = spherDynPVS.GetRadVel  ();
    }
    else
    {
      // Compute the GeoC Spherical Eq CoOrds:
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
             "%02.0lf:%02.0lf:%02.3lf\t"
           "%c%02.0lf %02.0lf %02.3lf\t"
             "%.2lf\t%.2lf\t%.3lf\n", 
           utc.m_year, utc.m_month, utc.m_day,
           utc.m_hour, utc.m_min,   utc.m_sec,
           get<0>(hms).Magnitude(), get<1>(hms).Magnitude(),
           get<2>(hms).Magnitude(),
           get<0>(dms) < 0  ? '-' : get<0>(dms) > 0 ? '+' : ' ',
           get<1>(dms).Magnitude(), get<2>(dms).Magnitude(),
           get<3>(dms).Magnitude(),
           alphaDotSD.Magnitude(),  deltaDotSD.Magnitude(),
           radVel.Magnitude());
  }
  // All Done!
  return 0;
}
