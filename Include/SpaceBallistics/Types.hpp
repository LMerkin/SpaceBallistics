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
    (Mass,  kg),
    (Len,   m,   (km, 1000.0)),
    (Time,  sec),

    // Angles: When expressed in "rad", angles are directly convertible to a
    // dimension-less "double";  other units must first be converted into "rad":
    (Angle,  rad,                       // Radians
      (deg,  Pi<double> / 180.0   ),    // Degrees
      (amin, Pi<double> / 10800.0 ),    // Arc-Minutes
      (asec, Pi<double> / 648000.0)     // Arc-Seconds
    )
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

  // Mass Rate (kg/sec):
  using MassRate   = decltype(1.0_kg   / 1.0_sec);

  // MoI Change Rate (kg*m^2 / sec):
  using MoIRate    = decltype(MoI(1.0) / 1.0_sec);
}
// End namespace SpaceBallistics


namespace DimTypes
{
  using namespace SpaceBallistics;

  //=========================================================================//
  // Conversion of Angles:                                                   //
  //=========================================================================//
  // Explicily allow conversion of "Angle[_rad]" into "double":
  template<>
  constexpr inline
  DimQ
  <
    Bits::DimExp(unsigned(DimsE::Angle)),
    Bits::MkUnit(unsigned(DimsE::Angle), 0),   // 0=rad: the main unit
    double
  >
  ::operator double() const
    { return Magnitude(); }
}
// End namespace DimTypes
