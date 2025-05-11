// vim:ts=2:et
//===========================================================================//
//                            "Tests/ERM_Test.cpp":                          //
//        Sun Longitude at Noon for a Year Containing a Leap Second          //
//===========================================================================//
#include "SpaceBallistics/CoOrds/TimeScales.hpp"
#include "SpaceBallistics/CoOrds/TopoCentricCOSes.h"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"
#include "SpaceBallistics/PhysEffects/DE440T.h"
#include "SpaceBallistics/PhysEffects/EarthRotationModel.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"

int main()
{
  using namespace SpaceBallistics;

  // Consts:
  constexpr int      Year   = 2015;   // With a Leap Second on June 30th
  constexpr Time_jyr EpochY (double(Year) + 0.5);
  constexpr TT       EpochTT(EpochY);

  // Construct the EarthRotationModel for that Year (centered at the middle
  // of the Year), using both the Analytical and the DE440T Nutations:
  //
  constexpr EarthRotationModel ERMAn(EpochY);
  EarthRotationModel           ERMDE(EpochTT);

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

    // ITRS position of the Sun, via the ERMs:
    PosKV_ITRS<Body::Sun>    itrsSunAn = ERMAn.ToITRS(tt, geoSun);
    PosKV_ITRS<Body::Sun>    itrsSunDE = ERMDE.ToITRS(tt, geoSun);

    // The TopoC (with ITRS Axes) positions of the Sun:
    PosKV<TCOS0, Body::Sun>  topSunAn  = ToTopoC<TCOS0>(itrsSunAn);
    PosKV<TCOS0, Body::Sun>  topSunDE  = ToTopoC<TCOS0>(itrsSunDE);

    // Convert the "topSun" co-ords into Spgerical ones: For the moment, we only
    // need the Lambda. This is the longitude of the Sun at Noon (should be aro-
    // und 0, since our location is in the XZ plane):
    SpherPV<TCOS0, Body::Sun> spherTopSunAn(topSunAn);
    SpherPV<TCOS0, Body::Sun> spherTopSunDE(topSunDE);

    // Output the Longitudes of the Sun, in Degrees:
    printf("%d\t%.9lf\t%.9lf\n",
           dc,
           To_Angle_deg(spherTopSunAn.GetAlpha()).Magnitude(),
           To_Angle_deg(spherTopSunDE.GetAlpha()).Magnitude());
  }
  // XXX: It appears that the discrepancy between the Sun Longitudes computed
  // with Analytical and DE440T ERMs is about 1e-6 deg, ie ~ 0".0036.
  // This is perfectly acceptable.
  return 0;
}
