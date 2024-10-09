// vim:ts=2:et
//===========================================================================//
//                        "SpaceBallistics/Types.hpp":                       //
//        Declarations of Principal Types used in "SpaceBallistics"          //
//===========================================================================//
#pragma  once
#include <DimTypes/DimTypes.hpp>
#include <SpaceBallistics/CoOrds/Vector3D.hpp>

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
    (Len,      m,   (km,     1000.0)),    // Length
    (Time,     sec, (day,   86400.0)),    // Time
    (AbsTemp,  K),                        // Absolute Temperature

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

  // Pressure:
  using Pressure = decltype(Force() / Area(1.0));

  // Angular Velocity and Angular Acceleration:
  using AngVel   = decltype(1.0    / 1.0_sec);
  using AngAcc   = decltype(1.0    / Sqr(1.0_sec));

  // Angular ("Kinetic") Momentum:
  using AngMom   = decltype(MoI()  * AngVel());

  // Rotational Moment of Force ("Torque"):
  using Torq     = decltype(Len()  * Force());

  // Mass Rate (kg/sec):
  using MassRate = decltype(Mass() / 1.0_sec);

  // MoI Change Rate (kg*m^2 / sec):
  using MoIRate  = decltype(MoI()  / 1.0_sec);

  // Gravitational Field Constant:
  using GM       = decltype(Cube(1.0_m) / Sqr(1.0_sec));

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

  //-------------------------------------------------------------------------//
  // Mechanical Vectors:                                                     //
  //-------------------------------------------------------------------------//
  // Macro for declaring a DimQ Vector (or a diagonal Tensor):
# ifdef DCL_VEC
# undef DCL_VEC
# endif
# define DCL_VEC(T) \
  template<typename COS> \
  using T##V = Vector3D<T, COS>;

  DCL_VEC(Len)      // Position Vector ("Radius-Vector")
  DCL_VEC(Vel)      // Velocity Vector
  DCL_VEC(Acc)      // Acceleration Vector
  DCL_VEC(Force)    // Force Vector
  DCL_VEC(AngVel)   // Angular Velocity Vector
  DCL_VEC(AngAcc)   // Angular Acceleration Vector
  DCL_VEC(AngMom)   // Angular ("Kinetic") Momentum Vector
  DCL_VEC(Torq)     // Rotational Moment of Force (Torque) Vector

  // Alias for the Position vector: "LenV" -> "PosV":
  template<typename COS>
  using PosV = LenV<COS>;

  // NB:  The MoI and its Rate of Change are in general not Vectors, but rather,
  // 3*3 Tensors. XXX: For the moment, we only consider those tensors in their
  // principal axes, so they have a diagonal form  and represented by Vectors:
  //
  DCL_VEC(MoI)      // Moments of Inertia
  DCL_VEC(MoIRate)  // MoI Change Rates
# undef DCL_VEC

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

