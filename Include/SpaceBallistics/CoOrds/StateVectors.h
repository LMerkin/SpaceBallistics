// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/StateVectors.h":                //
//                           Kinematic State Vectors                         //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // Fwd Declarations:                                                       //
  //=========================================================================//
  class BaryCentricCOS;
  class 

  //=========================================================================//
  // "StdKSV" Class:                                                         //
  //=========================================================================//
  // Kinematic State Vector for Steady Motion  (of some fixed embedded point or
  // the CoM, if the latter is considered to be fixed or slowly-moving), in the
  // specified COS. Can be used on its own if no Rotation Motion is considered,
  // or as part of a full 6D StateVector:
  //
  template<typename COS>
  struct StdKSV
  {
    using TimeScale = typename COS::TimeScale;

    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    TimeScale  const m_t;   // FIXME: May depend on the COS
    PosKV<COS> const m_r;  // (x,y,z) co-ords of the fixed point
    VelKV<COS> const m_v;  // (x_dot,  y_dot,  z_dot) velocity components

    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    constexpr StdKSV(TimeScale a_t, PosV<COS> const& a_r, VelV<COS> const& a_v)
    : m_t(a_t),
      m_r(a_r),
      m_v(a_v)
    {}
  };

  //=========================================================================//
  // "RotKSV" Struct:                                                        //
  //=========================================================================//
  // Kinematic State Vector for Rotational Motion (around some axis); the COS is
  // normally some "fixed" one, but it could also be a SnapShot of the Embedded
  // Co-Ords System (eg in the Control Theory problems).   XXX: Again, all data
  // flds are public, which is an easy-to-use but somewhat unsafe solution:
  //
  template<typename COS>
  struct RotKSV
  {
    using TimeScale = typename COS::TimeScale;

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Euler's Angles: (Pitch, Yaw, Roll). NB: They do not make a vector in the
    // mathematical sense; sometimes they are not well-defined  (eg for a snap-
    // shot of the EmbeddedCOS, the current Euler's angles would normally be 0):
    TimeScale    const m_t;
    Angle        const m_pitch;
    Angle        const m_yaw;
    Angle        const m_roll;

    // It is a good idea to memoise the Coses and Sines of Euler's Angles, for
    // use in further computations:
    double       const m_cosP;  // cos(Pitch)
    double       const m_sinP;  // sin(Pitch)
    double       const m_cosY;  // cos(Yaw)
    double       const m_sinY;  // sin(Yaw)
    double       const m_cosR;  // cos(Roll)
    double       const m_sinR;  // sin(Roll)

    // Euler's Angles Derivatives. In fact, they are uniquely determined by the
    // "m_omega" (below), and vice versa:
    AngVel       const m_pitchDot;
    AngVel       const m_yawDot;
    AngVel       const m_rollDot;

    // Vector of Angular Velocity. NB: This is NOT the vector of derivatives
    // (Pitch_dot, Yaw_dot, Roll_dot); the latter are given separately below:
    AngVelV<COS> const m_omega;

    //=======================================================================//
    // Non-Default Ctors:                                                    //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // With the Angular Velocity vector ("omega"):                           //
    //-----------------------------------------------------------------------//
    constexpr RotKSV
    (
      TimeScale           a_t,
      Angle               a_pitch,
      Angle               a_yaw,
      Angle               a_roll,
      AngVelV<COS> const& a_omega
    )
    : m_t       (a_t),
      m_pitch   (a_pitch),
      m_yaw     (a_yaw),
      m_roll    (a_roll),
      m_cosP    (Cos(m_pitch)),
      m_sinP    (Sin(m_pitch)),
      m_cosY    (Cos(m_yaw)),
      m_sinY    (Sin(m_yaw)),
      m_cosR    (Cos(m_roll)),
      m_sinR    (Sin(m_roll)),
      m_pitchDot(m_sinY * a_omega[0] - m_cosY * a_omega[1]),
      m_yawDot
      (
        // If "cosP" is close to 0, "yawDot" has a singularity; this normally
        // occurs during the vertical ascent, so we assume yawDot = 0 in this
        // case:
        (Abs(m_cosP) < Tol)
        ? AngVel(0.0)
        : a_omega[2] -
          (m_cosY * a_omega[0] + m_sinY * a_omega[1]) * m_sinP / m_cosP
      ),
      m_rollDot
      (
        // If "cosP" is close to 0, we assume that "rollDot" is +-omega[2], ie
        // omega[2] = omega_z is in this case translated into "rollDot",   not
        // into "yawDot" (whch remains 0, see above):
        (Abs(m_cosP) < Tol)
        ? (m_sinP > 0.0) ? a_omega[2] : (-a_omega[2])
        : (m_cosY * a_omega[0] + m_sinY * a_omega[1]) / m_cosP
      ),
      m_omega   (a_omega)
    {}

    //-----------------------------------------------------------------------//
    // With the Euler's Angles Derivatives:                                  //
    //-----------------------------------------------------------------------//
    constexpr RotKSV
    (
      TimeScale a_t,
      Angle     a_pitch,
      Angle     a_yaw,
      Angle     a_roll,
      AngVel    a_pitch_dot,
      AngVel    a_yaw_dot,
      AngVel    a_roll_dot
    )
    : m_t       (a_t),
      m_pitch   (a_pitch),
      m_yaw     (a_yaw),
      m_roll    (a_roll),
      m_cosP    (Cos(m_pitch)),
      m_sinP    (Sin(m_pitch)),
      m_cosY    (Cos(m_yaw)),
      m_sinY    (Sin(m_yaw)),
      m_cosR    (Cos(m_roll)),
      m_sinR    (Sin(m_roll)),
      m_pitchDot(a_pitch_dot),
      m_yawDot  (a_yaw_dot),
      m_rollDot (a_roll_dot),
      m_omega
      {
         m_sinY * m_pitchDot + m_cosP * m_cosY * m_rollDot,
        -m_cosY * m_pitchDot + m_cosP * m_sinY * m_rollDot,
         m_yawDot            + m_sinP          * m_rollDot
      }
    {}
  };

  //=========================================================================//
  // "KSV6D": Kinematic State Vector with 6 Degrees of Freedom:              //
  //=========================================================================//
  template<typename COS>
  class KSV6D: public StdKSV<COS>, public RotKSV<COS>
  {
  public:
    //=======================================================================//
    // Data Flds: Made directly visible:                                     //
    //=======================================================================//
    // XXX: Although both parent classes have the "m_t" fld, we do not use the
    // virtual inheritance here for the sake of simplicity and "constexpr":
    //
    using StdKSV<COS>::m_t;     // RotKSV<COS>::m_t will be the same
    using StdKSV<COS>::m_r;
    using StdKSV<COS>::m_v;
    using StdKSV<COS>::m_acc;

    using RotKSV<COS>::m_pitch;
    using RotKSV<COS>::m_yaw;
    using RotKSV<COS>::m_roll;
    using RotKSV<COS>::m_cosP;
    using RotKSV<COS>::m_sinP;
    using RotKSV<COS>::m_cosY;
    using RotKSV<COS>::m_sinY;
    using RotKSV<COS>::m_cosR;
    using RotKSV<COS>::m_sinR;
    using RotKSV<COS>::m_pitchDot;
    using RotKSV<COS>::m_yawDot;
    using RotKSV<COS>::m_rollDot;
    using RotKSV<COS>::m_omega;

    //=======================================================================//
    // Non-Default Ctors:                                                    //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // With the Angular Velocity vector ("omega"):                           //
    //-----------------------------------------------------------------------//
    constexpr KSV6D
    (
      TimeScale           a_t,
      PosV<COS>    const& a_r,
      VelV<COS>    const& a_v,
      Angle               a_pitch,
      Angle               a_yaw,
      Angle               a_roll,
      AngVelV<COS> const& a_omega
    )
    : StdKSV<COS>(a_t, a_r,     a_v),
      RotKSV<COS>(a_t, a_pitch, a_yaw, a_roll, a_omega)
    {}

    //-----------------------------------------------------------------------//
    // With the Euler's Angles Derivatives:                                  //
    //-----------------------------------------------------------------------//
    constexpr KSV6D
    (
      TimeScale           a_t,
      PosV<COS>    const& a_r,
      VelV<COS>    const& a_v,
      Angle               a_pitch,
      Angle               a_yaw,
      Angle               a_roll,
      AngVel              a_pitch_dot,
      AngVel              a_yaw_dot,
      AngVel              a_roll_dot
    )
    : StdKSV<COS>(a_t, a_r,     a_v),
      RotKSV<COS>(a_t, a_pitch, a_yaw, a_roll, a_pitch_dot, a_yaw_dot,
                  a_roll_dot)
    {}
  };
}
// End namespace SpaceBallistics
