// vim:ts=2:et
//===========================================================================//
//                            "PrecNutTest.cpp":                             //
//             Test for Long-Term Earth Precession and Nutations             //
//===========================================================================//
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"
#include "SpaceBallistics/CoOrds/SphericalPV.hpp"

int main()
{
  using namespace SpaceBallistics;

  // Full Precession Cycle from 2000.0 (the exact dates do not matter):
  constexpr Time_jyr From(2000.0);
  constexpr Time_jyr Step(1.0);
  constexpr Time_jyr To  = From + Time_jyr(26'000.0);

  // The Polar ITRS vector of unit length:
  constexpr PosKV<ITRS> TerrPolarVec { 0.0_km, 0.0_km, 1.0_km };

  for (Time_jyr t = From; t <= To; t += Step)
  {
    // ITRS orientation @ t:
    EarthRotationModel era{ t };

    // Convert the PolarVec to GCRS. For that, we have to convert "t" to TT,
    // although it does not matter here, since we only need  the Earth Axis
    // orientation:
    PosKV<GCRS> celestPolarVec = era.ToGCRS(TT{t}, TerrPolarVec);

    // Convert  the  "celestPolarVec" into Spherical Co-Ords:
    SphericalPV<GCRS> celestPolarSph(celestPolarVec, VelKV<GCRS>());

    // Output the "celestPolarVec" projections onto the GCRS Equatorial (XY)
    // plane, as well as the RA and Decl:
    printf("%.1lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n",
           t.Magnitude(),
           celestPolarVec.x()       .Magnitude(),
           celestPolarVec.y()       .Magnitude(),
           celestPolarSph.GetAlpha().Magnitude(),
           celestPolarSph.GetDelta().Magnitude());
  }
  return 0;
}
