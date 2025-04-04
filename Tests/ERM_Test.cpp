// vim:ts=2:et
//===========================================================================//
//                            "Tests/ERM_Test.cpp":                          //
//        Sun Longitude at Noon for a Year Containing a Leap Second          //
//===========================================================================//
#include "SpaceBallistics/CoOrds/TimeScales.hpp"
#include "SpaceBallistics/CoOrds/TopoCentricCOSes.h"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"
#include "SpaceBallistics/PhysEffects/DE440T.h"
#include "SpaceBallistics/PhysEffects/EarthRotationModel.h"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"

int main()
{
  using namespace SpaceBallistics;

  // Consts:
  constexpr int Year = 2015;    // With a Leap Second on June 30th

  // Construct the EarthRotationModel for that Year (centered at the middle
  // of the Year):
  EarthRotationModel ERM { Time_jyr(double(Year) + 0.5) };

  // TopoCRotCOS of our Observer's position (long=0, lat=0, h=0):
  constexpr static  Location<Body::Earth>   Observer(0.0_rad, 0.0_rad, 0.0_m);
  using     TCOS0 = TopoCRotCOS<Body::Earth, &Observer>;

  int dc = 0;
  for (int m = 1; m <= 12; ++m)
  for (int d = 1; d <= UTC::DaysInMonth(Year, m); ++d, ++dc)
  {
    // First, the UTC (@ Noon!)
    UTC utc(Year, m, d, 12, 0, 0.0);

    // Get the TT and TDB from the above UTC:
    TT  tt (utc);
    TDB tdb(tt);

    // Get the ICRS/BCRS (BaryC Equatorial J2000.0) Co-Ords of the Sun for this
    // time instant. The velocity is not required:
    PosKV_BCRS<Body::Sun>    posSun;
    DE440T::GetPlanetBarEqPV<Body::Sun>  (tdb, &posSun,   nullptr);

    // Same for Earth:
    PosKV_BCRS<Body::Earth>  posEarth;
    DE440T::GetPlanetBarEqPV<Body::Earth>(tdb, &posEarth, nullptr);

    // GeoC Equatorial position of the Sun:
    PosKV_GCRS<Body::Sun>    geoSun  = posSun - posEarth;

    // ITRS position of the Sun, via the ERM:
    PosKV_ITRS<Body::Sun>    itrsSun = ERM.ToITRS(tt, geoSun);

    // The TopoC (with ITRS Axes) position of the Sun:
    PosKV<TCOS0, Body::Sun>  topSun  = ToTopoC<TCOS0>(itrsSun);

    // Convert the "topSun" co-ords into Spgerical ones: For the moment, we only
    // need the Lambda. This is the longitude of the Sun at Noon (should be aro-
    // und 0, since our location is in the XZ plane):
    SpherPV<TCOS0, Body::Sun> spherTopSun(topSun);

    // Output the Longitude of the Sun, in Degrees:
    printf("%d\t%.6lf\n",
           dc, To_Angle_deg(spherTopSun.GetAlpha()).Magnitude());
  }
  return 0;
}
