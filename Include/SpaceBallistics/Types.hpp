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
    (Mass,     kg),
    (Len,      m,   (km, 1000.0)),
    (Time,     sec),

    // Angles: When expressed in "rad", angles are directly convertible to a
    // dimension-less "double";  other units must first be converted into "rad":
    (Angle,    rad,                       // Radians
      (deg,    Pi<double> / 180.0   ),    // Degrees
      (arcMin, Pi<double> / 10800.0 ),    // Arc-Minutes
      (arcSec, Pi<double> / 648000.0)     // Arc-Seconds
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

  // Moment of Inertia (kg*m^2):
  using MoI      = decltype(Mass() * Area(1.0));

  // Velocity (m/sec) and Acceleration (m/sec^2):
  using Vel      = decltype(Len()  / 1.0_sec);
  using Acc      = decltype(Len()  / Sqr(1.0_sec));

  // Force (N = kg*m/sec^2):
  using Force    = decltype(Mass() * Acc());

  // Angular Velocity and Angular Acceleration:
  using AngVel   = decltype(1.0    / 1.0_sec);
  using AngAcc   = decltype(1.0    / Sqr(1.0_sec));

  // Angular ("Kinetic") Momentum:
  using AngMom   = decltype(MoI()  * AngVel());

  // Rotational Moment of Force ("Torque"):
  using Torq     = decltype(Len()  * Force());

  // Standard Gravity (m/sec^2):
  constexpr inline Acc g0 = Acc(9.80665);

  // Mass Rate (kg/sec):
  using MassRate = decltype(Mass() / 1.0_sec);

  // MoI Change Rate (kg*m^2 / sec):
  using MoIRate  = decltype(MoI()  / 1.0_sec);

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
  // 3D Vectors, Parameterised by the CoOrd System (COS):                    //
  //=========================================================================//
  // XXX: They are just "std::array"s of size 3. All elements are initialised
  // to 0 by default (via the default "DimQ" Ctor).  The vectors are mutable,
  // which may be somewhat unsafe(?).
  // Macro for declaring a Vector (or a diagonal Tensor):
# ifdef DCL_VEC
# undef DCL_VEC
# endif
  // The following macro creates the {T}{Sfx} Vec3 templated by COS, for the
  // type "T"; the "Sfx" is "V" for a Vector and "T" for a Tensor:
  //
# define DCL_VEC(T, Sfx) \
  template<typename COS> \
  using T##Sfx = std::array<T, 3>;

  DCL_VEC(Len,     V)  // Position Vector ("Radius-Vector")
  DCL_VEC(Vel,     V)  // Velocity Vector
  DCL_VEC(Acc,     V)  // Acceleration Vector
  DCL_VEC(Force,   V)  // Force Vector
  DCL_VEC(AngVel,  V)  // Angular Velocity Vector
  DCL_VEC(AngAcc,  V)  // Angular Acceleration Vector
  DCL_VEC(AngMom,  V)  // Angular ("Kinetic") Momentum Vector
  DCL_VEC(Torq,    V)  // Rotational Moment of Force (Torque) Vector

  // Alias for the Position vector: "LenV" -> "PosV":
  template<typename COS>
  using PosV = LenV<COS>;

  // NB:  The MoI and its Rate of Change are in general not Vectors, but rather,
  // 3*3 Tensors. XXX: For the moment, we only consider those tensors in their
  // principal axes, so they have a diagonal form. Hence the suffix is "T", not
  // "V", although the internal rep is still "Vec3":
  //
  DCL_VEC(MoI,     T)   // Moments of Inertia
  DCL_VEC(MoIRate, T)   // MoI Change Rates
# undef DCL_VEC

  //=========================================================================//
  // Computation Tolerances:                                                 //
  //========================================================================//
  // (Searching for a better place to define them...):
  constexpr inline double Tol     = 100.0 * Eps<double>;
  constexpr inline double TolFact = 1.0   + Tol;
}
// End namespace SpaceBallistics
