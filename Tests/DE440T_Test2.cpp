// vim:ts=2:et
//===========================================================================//
//                        "Tests/DE440T_Test2.cpp":                          //
//===========================================================================//
#include <cstdio>
#include <cstdlib>

// Asserts must be enabled for this test:
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"
#include "SpaceBallistics/Utils.hpp"

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

    PosKVBarEq posE, posS;
    VelKVBarEq velE, velS;
    DE440T::GetPlanetBarEqPV<Body::Earth>(tdb, &posE, &velE);
    DE440T::GetPlanetBarEqPV<Body::Sun>  (tdb, &posS, &velS);

    // Compute the GeoEq PV vectors of the Sun:
    PosKVGeoEqFix posES
      (posS.x() - posE.x(), posS.y() - posE.y(), posS.z() - posE.z());
    VelKVGeoEqFix velES
      (velS.x() - velE.x(), velS.y() - velE.y(), velS.z() - velE.z());

    // Compute the GeoCentric Spherical Eq CoOrds:
    GeoCentricEqSpherPV spherPVS(posES, velES);

    auto hms = ToHMS(spherPVS.GetRA  ());
    auto dms = ToDMS(spherPVS.GetDecl());

    // Output:
    printf("%.4lf\t"
             "%02.0lf:%02.0lf:%02.3lf\t"
           "%c%02.0lf %02.0lf %02.3lf\n",
           tt.GetJD() .Magnitude(),
           get<0>(hms).Magnitude(), get<1>(hms).Magnitude(),
           get<2>(hms).Magnitude(),
           get<0>(dms) < 0  ? '-' : get<0>(dms) > 0 ? '+' : ' ',
           get<1>(dms).Magnitude(), get<2>(dms).Magnitude(),
           get<3>(dms).Magnitude());
  }
  // All Done!
  return 0;
}
