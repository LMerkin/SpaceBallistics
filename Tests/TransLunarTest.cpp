// vim:ts=2:et
//============================================================================//
//                         "Tests/TransLunarTest.cpp":                        //
//============================================================================//
#include "SpaceBallistics/CoOrds/EllipticOrbit.hpp"

using namespace SpaceBallistics;

namespace
{
  //==========================================================================//
  // Instantaneous Keplerian Lunar Orbit (in GCRS):                           //
  //==========================================================================//
  EllipticOrbit<GCRS, Body::Moon> GetInstMoonOrbit(TDB a_tdb)
  {
    PosKV_GCRS<Body::Moon> posM;
    VelKV_GCRS<Body::Moon> velM;
    DE440T::GetMoonGEqPV(a_tdb, &posM, &velM);

    return EllipticOrbit<GCRS,  Body::Moon>(a_tdb, posM, velM);
  }
}
 
//============================================================================//
// "main":                                                                    //
//============================================================================//
int main()
{
  constexpr Time_day From = Epoch_J2000;
  constexpr Time_day To   = From + To_Time_day(40.0_jyr);
  constexpr Time_day Step = 30.0_day;

  for (Time_day t = From; t <= To; t += Step)
  {
    TDB tdb{TT{t}};
    EllipticOrbit<GCRS, Body::Moon> orbitM = GetInstMoonOrbit(tdb);

    Angle_deg I     = To_Angle_deg(orbitM.I());
    Angle_deg Omega = To_Angle_deg(orbitM.Omega());
    Angle_deg omega = To_Angle_deg(orbitM.omega());
    Angle_deg M0    = To_Angle_deg(orbitM.M0());
    Time_day  T     = orbitM.T().GetJD();

    printf("%.3lf\t%.6lf\t%.6lf\t%.3lf\t%.9lf\t%.6lf\t%.6lf\t%.6lf\n",
           t         .Magnitude(),
           I         .Magnitude(),
           Omega     .Magnitude(),
           orbitM.a().Magnitude(),
           orbitM.e(),
           omega     .Magnitude(),
           M0        .Magnitude(),
           T         .Magnitude());
  }
  return 0;
}
