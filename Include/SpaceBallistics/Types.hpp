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
    double,    8,                         // Up to 8 dims
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
  using Acc      = decltype(Vel()  / 1.0_sec);

  // Using AstroDynamical Units (km, km/sec):
  using LenK     = Len_km;
  using VelK     = decltype(LenK() / 1.0_sec);
  using AccK     = decltype(VelK() / 1.0_sec);

  // Force (N = kg*m/sec^2) and the resp "ForceK":
  using Force    = decltype(Mass() * Acc ());
  using ForceK   = decltype(Mass() * AccK());

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
  constexpr inline Acc      g0     = Acc(9.80665);
  constexpr inline AccK     g0K    = To_Len_km(g0);

  // "Pi"-Related Consts lifted to "Angle"s:
  constexpr inline Angle    PI     = Angle(Pi   <double>);
  constexpr inline Angle    TWO_PI = Angle(TwoPi<double>);
  constexpr inline Angle    PI_2   = Angle(Pi_2 <double>);

  //=========================================================================//
  // Helpers for "Angle"s and "DimLess" qtys:                                //
  //=========================================================================//
  // Explicit conversion of "Angle[_rad]" into "DimLess".
  // XXX: Using "Magnitude()" for efficiency, instead of division by "1.0_rad":
  //
  constexpr  DimLess To_DimLess(Angle a_angle)
    { return DimLess(a_angle.Magnitude());   }

  // Trignomometric Functions on "Angle"s. XXX: HOWEVER, they will hide similar
  // functions on "double"s -- but that may actually be a good idea:
  //
  constexpr double     Sin (Angle a_x)
    { return DimTypes::Sin (double(To_DimLess(a_x))); }

  constexpr double     Cos (Angle a_x)
    { return DimTypes::Cos (double(To_DimLess(a_x))); }

  constexpr double     Tan (Angle a_x)
    { return DimTypes::Tan (double(To_DimLess(a_x))); }

  // Elementary Unary functions of Non-Angular "DimLess" qtys:
  constexpr double ASin (DimLess a_x) { return DimTypes::ASin (double(a_x)); }
  constexpr double ACos (DimLess a_x) { return DimTypes::ACos (double(a_x)); }
  constexpr double ATan (DimLess a_x) { return DimTypes::ATan (double(a_x)); }
  constexpr double Exp  (DimLess a_x) { return DimTypes::Exp  (double(a_x)); }
  constexpr double Log  (DimLess a_x) { return DimTypes::Log  (double(a_x)); }
  constexpr double SinH (DimLess a_x) { return DimTypes::SinH (double(a_x)); }
  constexpr double CosH (DimLess a_x) { return DimTypes::CosH (double(a_x)); }
  constexpr double TanH (DimLess a_x) { return DimTypes::TanH (double(a_x)); }
  constexpr double ASinH(DimLess a_x) { return DimTypes::ASinH(double(a_x)); }
  constexpr double ACosH(DimLess a_x) { return DimTypes::ACosH(double(a_x)); }
  constexpr double ATanH(DimLess a_x) { return DimTypes::ATanH(double(a_x)); }
}
// End namespace SpaceBallistics

namespace DimTypes
{
  namespace SB = SpaceBallistics;

  //-------------------------------------------------------------------------//
  // Explicit conversions of "Angle[_rad]" into "double":                    //
  //-------------------------------------------------------------------------//
  // XXX: This function is in the form of "operator double", and therefore, it
  // must be defined in the same namespace ("DimTypes") as the original "DimQ"s:
  //
  template<>
  constexpr inline
  DimQ
  <
    SB::DimQ_Encs::DimExp(unsigned(SB::DimsE::Angle)),
    SB::DimQ_Encs::MkUnit(unsigned(SB::DimsE::Angle), 0), // 0=rad: fund unit
    SB::DimQ_RepT,
    SB::DimQ_MaxDims
  >
  ::operator double() const 
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle_rad>);
    return double(SB::To_DimLess(*this));
  }

  //-------------------------------------------------------------------------//
  // Similarly, converting "Angle_deg" into "double", via "Angle" (rad):     //
  //-------------------------------------------------------------------------//
  template<>
  constexpr inline
  DimQ
  <
    SB::DimQ_Encs::DimExp(unsigned(SB::DimsE::Angle)),
    SB::DimQ_Encs::MkUnit(unsigned(SB::DimsE::Angle), 1), // 1=deg
    SB::DimQ_RepT,
    SB::DimQ_MaxDims
  >
  ::operator double() const
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle_deg>);
    return double(SB::To_Angle(*this));
  }
}
// End namespace DimTypes

//===========================================================================//
// Support for using "std::format" with our "DimQ"s:                         //
//===========================================================================//
MK_DIMS_FMT(SpaceBallistics)

