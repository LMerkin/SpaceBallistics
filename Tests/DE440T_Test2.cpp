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
#include "SpaceBallistics/CoOrds/SphericalPV.hpp"
#include "SpaceBallistics/Utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cassert>

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  if (argc != 2 && argc != 3 && argc != 4)
  {
    cerr << "PARAMS: TT_From [TT_To [TimeStep_days]]" << endl;
    return 1;
  }
  // Get the JD_TT:
  TT   from    {Time_day(atof(argv[1]))};
  TT   to    = (argc >= 3) ?      TT{Time_day(atof(argv[2]))} : from;
  Time step  = To_Time((argc >= 4) ? Time_day(atof(argv[3]))  : 1.0_day);

  if (TDB(from) < DE440T::Bits::From || TDB(to) > DE440T::Bits::To ||
      from      > to                 || !IsPos(step))
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
    // FIXME: Replace the following with type-safe vector subtractions!
    PosKV_GCRS<Body::Sun> posES = posS - posE;
    VelKV_GCRS<Body::Sun> velES = velS - velE;

    // Compute the GeoCentric Spherical Eq CoOrds:
    GeoCentricEqSpherPV<Body::Sun> spherPVS(posES, velES);

    auto hms = ToHMS(spherPVS.GetAlpha());
    auto dms = ToDMS(spherPVS.GetDelta());

    // Compute the RA and Decl "Dots" (in Seconds / ArcSeconds per Day):
    auto alphaDot = To_Time_day(To_Angle_ss    (spherPVS.GetAlphaDot()));
    auto deltaDot = To_Time_day(To_Angle_arcSec(spherPVS.GetDeltaDot()));

    // Output:
    printf("%.4lf\t"
             "%02.0lf:%02.0lf:%02.3lf\t"
           "%c%02.0lf %02.0lf %02.3lf\t"
             "%.2lf\t%.2lf\t%.3lf\n",
           tt.GetJD() .Magnitude(),
           get<0>(hms).Magnitude(), get<1>(hms).Magnitude(),
           get<2>(hms).Magnitude(),
           get<0>(dms) < 0  ? '-' : get<0>(dms) > 0 ? '+' : ' ',
           get<1>(dms).Magnitude(), get<2>(dms).Magnitude(),
           get<3>(dms).Magnitude(),
           alphaDot.Magnitude(),    deltaDot.Magnitude(),
           spherPVS.GetRadVel().Magnitude());

    // XXX: Although we explicitly reset the NDEBUG flag at the beginning, in
    // CLang it may still be set (???), so use the following guard to prevent
    // compile-time warnings.   However, this will prevent the following test
    // from being done:
#   ifndef NDEBUG
    // Get back to the GeoEq PV vectors:
    auto  [posES1, velES1] = spherPVS.GetPVVectors();
    assert(posES1.x().ApproxEquals(posES.x()) &&
           posES1.y().ApproxEquals(posES.y()) &&
           posES1.z().ApproxEquals(posES.z()) &&
           velES1.x().ApproxEquals(velES.x()) &&
           velES1.y().ApproxEquals(velES.y()) &&
           velES1.z().ApproxEquals(velES.z()));
#   endif
  }
  // All Done!
  return 0;
}
