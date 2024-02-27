// vim:ts=2:et
//===========================================================================//
//                                "Types.hpp":                               //
//        Declarations of Principal Types used in "SpaceBallistics"          //
//===========================================================================//
#pragma  once
#include <DimTypes/DimTypes.hpp>

namespace SpaceBallistics
{
  using namespace DimTypes;

  //=========================================================================//
  // Basic Dimension Types and Units:                                        //
  //=========================================================================//
  DECLARE_DIMS
  (
    double,
    (Mass, kg),
    (Len,  m, (km, 1000.0)),
    (Time, sec)
  )

  //=========================================================================//
  // Derived Dimension Types and Values:                                     //
  //=========================================================================//
  // Area (m^2) and Volume (m^3):
  using Area       = decltype(Sqr (1.0_m));
  using Vol        = decltype(Cube(1.0_m));

  // Density (ie Volume Density, kg/m^3) and Surface (Area) Density (kg/m^2):
  using Density    = decltype(1.0_kg / Vol (1.0));
  using SurfDens   = decltype(1.0_kg / Area(1.0));

  // Moment of Inertia (kg*m^2):
  using MoI        = decltype(1.0_kg * Area(1.0));

  // Velocity (m/sec) and Acceleration (m/sec^2):
  using Vel        = decltype(1.0_m  / 1.0_sec);
  using Acc        = decltype(1.0_m  / Sqr(1.0_sec));

  // Angular Velocity and Angular Acceleration:
  using AngVel     = decltype(1.0    / 1.0_sec);
  using AngAcc     = decltype(1.0    / Sqr(1.0_sec));

  // Standard Gravity (m/sec^2):
  constexpr inline Acc g0 = Acc(9.80665);

  // Force (N = kg*m/sec^2):
  using Force      = decltype(1.0_kg * Acc(1.0));
}
// End namespace SpaceBallistics
