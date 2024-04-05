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
  using Area       = decltype(Sqr (Len()));
  using Vol        = decltype(Cube(Len()));

  // Density (ie Volume Density, kg/m^3) and Surface (Area) Density (kg/m^2):
  using Density    = decltype(Mass() / Vol (1.0));
  using SurfDens   = decltype(Mass() / Area(1.0));

  // Moment of Inertia (kg*m^2):
  using MoI        = decltype(Mass() * Area(1.0));

  // Velocity (m/sec) and Acceleration (m/sec^2):
  using Vel        = decltype(Len()  / 1.0_sec);
  using Acc        = decltype(Len()  / Sqr(1.0_sec));

  // Angular Velocity and Angular Acceleration:
  using AngVel     = decltype(1.0    / 1.0_sec);
  using AngAcc     = decltype(1.0    / Sqr(1.0_sec));

  // Standard Gravity (m/sec^2):
  constexpr inline Acc g0 = Acc(9.80665);

  // Force (N = kg*m/sec^2):
  using Force      = decltype(Mass() * Acc());

  // Mass Rate (kg/sec):
  using MassRate   = decltype(Mass() / 1.0_sec);

  // MoI Change Rate (kg*m^2 / sec):
  using MoIRate    = decltype(MoI()  / 1.0_sec);

  //-------------------------------------------------------------------------//
  // Powers of "Len" and their Time Derivatives: Widely used:                //
  //-------------------------------------------------------------------------//
  using Len2       = Area;
  using Len3       = Vol;
  using Len4       = decltype(Sqr(Area()));
  using Len5       = decltype(Len4() * Len());
  using Len6       = decltype(Sqr(Vol ()));

  using Len2Rate   = decltype(Len2() / 1.0_sec);
  using Len3Rate   = decltype(Len3() / 1.0_sec);
  using VolRate    = Len3Rate;
  using Len4Rate   = decltype(Len4() / 1.0_sec);
  using Len5Rate   = decltype(Len5() / 1.0_sec);

  //=========================================================================//
  // 3D Vectors, Parameterised by the CoOrd System (COS):                    //
  //=========================================================================//
}
// End namespace SpaceBallistics


namespace DimTypes
{
  namespace SB = SpaceBallistics;

  //=========================================================================//
  // Conversion of Angles:                                                   //
  //=========================================================================//
  // Explicily allow conversion of "Angle[_rad]" into "double":
  template<>
  constexpr inline
  DimQ
  <
    Bits::DimExp(unsigned(SB::DimsE::Angle)),
    Bits::MkUnit(unsigned(SB::DimsE::Angle), 0),   // 0=rad: the main unit
    double
  >
  ::operator double() const
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle>   &&
                  std::is_same_v<ThisDimQ, SB::Angle_rad>);
    return Magnitude();
  }

  // Similarly, converting "Angle_deg" into "double", via "Angle_rad":
  template<>
  constexpr inline
  DimQ
  <
    Bits::DimExp(unsigned(SB::DimsE::Angle)),
    Bits::MkUnit(unsigned(SB::DimsE::Angle), 1),   // 1=deg
    double
  >
  ::operator double() const
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle_deg>);
    return double(SB::To_Angle_rad(*this));
  }
}
// End namespace DimTypes
