// vim:ts=2:et
//===========================================================================//
//                        "SpaceBallistics/Types.hpp":                       //
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
    double,,
    (Mass,     kg),                       // Mass
    (Len,      m,     (km,  1000.0)),     // Length

    // Time:
    (Time,     sec,
      (day,    86400.0),
      (tyr,    31'556'925.9747),          // Tropical Year 1900 (CGPM 1960) ~=
                                          //   365.2422_day
      (jyr,    31'557'600.0)              // Julian   Year = 365.25_day
    ),

    (AbsTemp,  K),                        // Absolute Temperature

    // Angles: When expressed in "rad", angles are directly convertible to a
    // dimension-less "double";  other units must first be converted into "rad":
    (Angle,    rad,                       // Radians
      // (deg, min, sec):
      (deg,    Pi<double> / 180.0   ),    // Degrees
      (arcMin, Pi<double> / 10800.0 ),    // Arc-Minutes
      (arcSec, Pi<double> / 648000.0),    // Arc-Seconds
      // (hh,  mm,  ss):
      (hh,     Pi<double> / 12.0    ),    // Time-Hours   as Angles
      (mm,     Pi<double> / 720.0   ),    // Time-Minutes as Angles
      (ss,     Pi<double> / 43200.0 ),    // Time-Seconds as Angles
    )
  )

  //=========================================================================//
  // Derived Dimension Types and Values:                                     //
  //=========================================================================//
  // Area (m^2) and Volume (m^3):
  using Area     = decltype(Sqr (Len()));
  using Vol      = decltype(Cube(Len()));

  // Density (ie Volume Density, kg/m^3) and Surface (Area) Density (kg/m^2):
  using Density  = decltype(Mass() / Vol (1.0));
  using SurfDens = decltype(Mass() / Area(1.0));

  // For compressible gases: Rate of Density change:
  using DensRate = decltype(Density() / 1.0_sec);

  // Moment of Inertia (kg*m^2):
  using MoI      = decltype(Mass() * Area(1.0));

  // Velocity (m/sec) and Acceleration (m/sec^2):
  using Vel      = decltype(Len()  / 1.0_sec);
  using Acc      = decltype(Len()  / Sqr(1.0_sec));

  // Using AstroDynamical Units (km, km/sec):
  using LenK     = Len_km;
  using VelK     = decltype(LenK() / 1.0_sec);

  // Force (N = kg*m/sec^2):
  using Force    = decltype(Mass() * Acc());

  // Pressure:
  using Pressure = decltype(Force() / Area(1.0));

  // Angular Velocity and Angular Acceleration:
  using AngVel   = decltype(Angle() / 1.0_sec);
  using AngAcc   = decltype(Angle() / Sqr(1.0_sec));

  // Angular ("Kinetic") Momentum:
  using AngMom   = decltype(MoI()  * AngVel());

  // Rotational Moment of Force ("Torque"):
  using Torq     = decltype(Len()  * Force());

  // Mass Rate (kg/sec):
  using MassRate = decltype(Mass() / 1.0_sec);

  // MoI Change Rate (kg*m^2 / sec):
  using MoIRate  = decltype(MoI()  / 1.0_sec);

  // Gravitational Field Constant (NB: km^3/sec^2):
  using GMK      = decltype(Cube(1.0_km) / Sqr(1.0_sec));

  // Energy:
  using Energy   = decltype(Force(1.0) * 1.0_m);

  //-------------------------------------------------------------------------//
  // Powers of "Len" and their Time Derivatives: Widely used:                //
  //-------------------------------------------------------------------------//
  using Len2     = Area;
  using Len3     = Vol;
  using Len4     = decltype(Sqr(Area()));
  using Len5     = decltype(Len4() * Len());
  using Len6     = decltype(Sqr(Vol ()));

  using Len2Rate = decltype(Len2() / 1.0_sec);
  using Len3Rate = decltype(Len3() / 1.0_sec);
  using VolRate  = Len3Rate;
  using Len4Rate = decltype(Len4() / 1.0_sec);
  using Len5Rate = decltype(Len5() / 1.0_sec);

  //=========================================================================//
  // Some Standard Constants:                                                //
  //=========================================================================//
  // Standard Gravity (m/sec^2):
  constexpr inline Acc      g0   = Acc(9.80665);

  // Standard Atnospheric Pressure at Sea Level:
  constexpr inline Pressure p0   = Pressure(101325.0);

  // Specific Gas Constant for Dry Air (J/K):
  constexpr inline auto     Rair = 287.0528 * Energy(1.0) / 1.0_K;

  //=========================================================================//
  // Computation Tolerances:                                                 //
  //=========================================================================//
  // (Searching for a better place to define them...):
  constexpr inline double Tol     = 100.0 * Eps<double>;
  constexpr inline double TolFact = 1.0   + Tol;
}
// End namespace SpaceBallistics

//===========================================================================//
// Support for using "std::format" with our "DimQ"s:                         //
//===========================================================================//
MK_DIMS_FMT(SpaceBallistics)

