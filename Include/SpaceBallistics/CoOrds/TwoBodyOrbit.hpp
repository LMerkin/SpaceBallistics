// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/TwoBodyOrbit.hpp":                 //
//              Kinematics of the Keplerian (Two-Body) Motion                //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "TwoBodyOrbit" Class:                                                   //
  //=========================================================================//
  template<typename COS, Body B = Body::UNDEFINED>
  class TwoBodyOrbit
  {
  public:
    //=======================================================================//
    // Orbit Types:                                                          //
    //=======================================================================//
    enum class Kind
    {
      Circular       = 0,
      Elliptic       = 1,
      Parabolic      = 2,
      Hyperbolic     = 3,
      RectElliptic   = 4,   // RectiLinear Ellipse
      RectParabolic  = 5,   // RectiLinear Parabola
      RectHyperbolic = 6    // RectiLinear Hyperbola
    };

  private:
    //=======================================================================//
    // Elliptic Keplerian Elements:                                          //
    //=======================================================================//
    Kind   m_kind;
    LenK   m_a;      // Semi-Major Axis
    double m_e;      // Eccentricity
    Angle  m_I;      // Orbital Plane Inclination (to GCRS XY=Eq Plane)
    Angle  m_Omega;  // Ascending Node Longitude  (from X)
    Angle  m_omega;  // Argument  of the PeriFocus
    TDB    m_T;      // PeriFocus Passage Time (in TDB for DE440T)

    // NB: From the 2-Body Problem, the Gravitational Fld Constant must be the
    // sum of such Constants for the Focal Body and for "B":
    constexpr static GMK Mu = DE440T::K<COS::BaseBody> + DE440T::K<B>;

    // Some Derived Elemennts:
    LenK   m_q;      // PeriFocus Distance
    LenK   m_Q;      // ApoFocus  Distance
    LenK   m_p;      // Focal Parameter
    decltype(Sqr(1.0_km) / 1.0_sec)
           m_h;      // Angular Momentum Constant
    double m_S;      // SqRt((1+e)/(1-e))
    AngVel m_n;      // Mean Motion (rad/sec !!!)
    Time   m_P;      // Orbital Period  (sec !!!)
    TDB    m_t0;     // Time of these Elements
    Angle  m_E0;     // Eccentric Anomaly (or its variants) @ m_t0
    Angle  m_M0;     // Mean Anomaly @ m_t0
    Angle  m_l0;     // Longitude from the Ascending Node   @ m_t0
    Angle  m_f0;     // True Anomaly @ m_t0

    // "TM": Matrix for Transforming "In-Orbit" Co-Ords [X', Y', Z'] into the
    // "COS", where:
    // X' axis: from the Focus to the PeriFocus;
    // Z' axis: from the Focus, normal to the orbital plane, so that the orbital
    //          motion is Counter-Clock-Wise when viewed from the end of Z';
    // Y' axis: such that [X',Y',Z'] is the right-oriented frame:
    //
    Mtx33  m_TM;     // Rotation Matrix: R = TM * R'
    Mtx33  m_invTM;  // Inverse (= Transposed)  "TM"

    //=======================================================================//
    // Default Ctor is for the internal use only:                            //
    //=======================================================================//
    TwoBodyOrbit() = default;

  public:
    //=======================================================================//
    // Non-Default Ctor: From Pos and Vel:                                   //
    //=======================================================================//
    TwoBodyOrbit
    (
      TDB                           a_t0,
      Vector3D<LenK, COS, B> const& a_pos,
      Vector3D<VelK, COS, B> const& a_vel
    )
    {
      //---------------------------------------------------------------------//
      // Radius-Vector and Velocity:                                         //
      //---------------------------------------------------------------------//
      // Radius-Vector must be > 0, whereas the Velocity may be 0 in case of a
      // RectElliptic orbit:
      LenK r = LenK(a_pos);
      assert(IsPos(r));
      VelK V = VelK(a_vel);

      //---------------------------------------------------------------------//
      // The Angular Momentum Vector:                                        //
      //---------------------------------------------------------------------//
      auto angMom = CrossProd(a_pos, a_vel);
      m_h         = angMom.EuclidNorm();

      // Degenetate Case: RectiLinear orbit:
      bool isRect = false;
      if (double(m_h / (r * V)) < DefaultTol<double>)
      {
        isRect = true;
        m_h    = decltype(m_h)(0.0);
      }

      //---------------------------------------------------------------------//
      // Inclination and the Longitude of the Ascending Node:                //
      //---------------------------------------------------------------------//
      bool Iis0  = false;
      bool IisPi = false;

      if (LIKELY(!isRect))
      {
        //-------------------------------------------------------------------//
        // Generic Case: A Planar Orbit:                                     //
        //-------------------------------------------------------------------//
        //  "k" is the unit vector orthogonal to the Orbital Plane. It is well-
        //  defined for non-RectLin orbits:
        Vector3D<DimLess, COS, B> k = angMom / m_h;

        // Orbit Inclination:
        m_I   = Angle(ACos(double(k.z())));
        assert(0.0_rad <= m_I && m_I <= PI);

        // Are we close to I=0 or I=Pi?
        constexpr DimLess Near1(1.0 - 1e-6);
        Iis0  = k.z() >  Near1;
        IisPi = k.z() < -Near1;

        if (UNLIKELY(Iis0))
          m_I = 0.0_rad;
        else
        if (UNLIKELY(IisPi))
          m_I = PI;

        // Ascending Node Longitude:
        // k.x = sin(I) * sin(Omega), k.y = -sin(I) * cos(Omega), sin(I) >= 0;
        // XXX:  For I=0 and I=Pi,  the Orbit is in the XY plane of the "COS",
        // "m_Omega" is undefined and is set to 0:
        m_Omega =
          (LIKELY(!(Iis0 || IisPi)))
          ? To2Pi(Angle(ATan2(k.x(), -k.y())))
          : 0.0_rad;

        //-------------------------------------------------------------------//
        // Longitude in orbit from the Ascending Node @ "m_t0":              //
        //-------------------------------------------------------------------//
        if (LIKELY(!(Iis0 || IisPi)))
        {
          double cosl0 =
            double((a_pos.x() * Cos(m_Omega)  +
                    a_pos.y() * Sin(m_Omega)) / r);
          m_l0 = Angle(ACos(cosl0));      // In [0 .. Pi] initially

          // Possibly adjust "l0" in the generic case:
          if (IsNeg(a_pos.z()))
            // Then, obviously, Pi < l0 < 2*Pi, so we must update "l0" but pre-
            // serve the "cosl0":
            m_l0 = TWO_PI - m_l0;
          assert(0.0_rad <= m_l0 && m_l0 < TWO_PI);
        }
        else
          // "m_l0" is undefined for I==0 and I==Pi; the motion occirs in the XY
          // plane of the "COS", so "m_l0" is set to the longitude from X, not
          // from the (non-existent) Ascending Node.
          // NB: This way, "m_l0" is defined even for Curcular orbits for which
          // the is no PeriFocus!
          m_l0 = To2Pi(Angle(ATan2(a_pos.y(), a_pos.x())));
      }
      else
      {
        //-------------------------------------------------------------------//
        // Degenerate Case: RectiLinear Orbit:                               //
        //-------------------------------------------------------------------//
        // For a RectLin orbit, "I" and "Omega" can also be defined, as Euler's
        // Angles of the radius-vector. But in this case, we may get I < 0:
        m_I     = Angle(ASin(double(a_pos.z() / r)));
        assert(-PI_2 <= m_I && m_I <= PI/2);
        m_Omega = Angle(ATan2(a_pos.y(), a_pos.x()));

        // In this case, "m_l0" is undefined, because there is no Ascending No-
        // de actually;  set it to 0:
        m_l0    = 0.0_rad;
      }

      //---------------------------------------------------------------------//
      // Semi-Major Axis, from the Energy Integral:                          //
      //---------------------------------------------------------------------//
      // V^2  = Mu * (2/r - 1/a) :
      auto invA = (2.0 / r - Sqr(V) / Mu);

      if (Abs(invA) < 1e-6 / r)     // XXX: Probably reasonable...
        invA = decltype(invA)(0.0); // Will be a Parabolic orbit

      // We can now provisionally set the Orbit Kind (except for the Circular
      // orbit which is as yet classified as Elliptic):
      m_kind  =
        IsPos(invA)
        ? (UNLIKELY(isRect) ? Kind::RectElliptic   : Kind::Elliptic)  :
        IsZero(invA)
        ? (UNLIKELY(isRect) ? Kind::RectParabolic  : Kind::Parabolic) :
          (UNLIKELY(isRect) ? Kind::RectHyperbolic : Kind::Hyperbolic);

      // Infinity for Parabolic, that's OK; for Hyperbolic and RecHyperbolic,
      // we still have m_a > 0:
      m_a = 1.0 / Abs(invA);
      assert(IsPos(m_a));

      //---------------------------------------------------------------------//
      // Focal Param:                                                        //
      //---------------------------------------------------------------------//
      // The following works in all cases; we just need a fix for possible
      // rounding errors in the Elliptic case:
      m_p = Sqr(m_h) / Mu;
      assert(IsZero(m_p) == isRect);

      if (m_kind == Kind::Elliptic)
      {
        assert(IsPos(m_p));
        // We must have (m_p <= m_a), so enforce it:
        m_p = std::min(m_p, m_a);
      }

      //---------------------------------------------------------------------//
      // The Eccentricity:                                                   //
      //---------------------------------------------------------------------//
      m_e = SqRt(1.0 - double(m_p / m_a));

      assert(0.0 <= m_e);
      assert((m_kind == Kind::Elliptic)     == (m_e < 1.0));
      assert((m_kind == Kind::Parabolic     || m_kind == Kind::RectElliptic ||
              m_kind == Kind::RectParabolic || m_kind == Kind::RectHyperbolic)
             == (m_e == 1.0));
      assert((m_kind == Kind::Hyperbolic)   == (m_e > 1.0));

      // The special case of the Circular orbit:
      if (m_e < 1e-6)
      {
        m_e    = 0.0;
        m_kind = Kind::Circular;
      }

      //---------------------------------------------------------------------//
      // Derived Elements:                                                   //
      //---------------------------------------------------------------------//
      // PeriFocus Distance: Always well-defined. Since we must always have
      // m_q <= r, enforce it:
      m_q = std::min(m_a * Abs(1.0 - m_e), r);

      // ApoFocus Distance: Explicitly set to +oo for non-periodic orbits:
      m_Q =
        (m_kind == Kind::Parabolic  || m_kind == Kind::RectParabolic ||
         m_kind == Kind::Hyperbolic || m_kind == Kind::RectHyperbolic)
        ? LenK(+Inf<double>)
        : std::max(m_a * (1.0 + m_e), r);     // Enforce Q >= r

      // A constant used in Elliptic and Hyperbolic variants of Kepler's Equat-
      // ion:
      m_S =
        (m_kind == Kind::Circular   || m_kind == Kind::Elliptic)
        ? SqRt ((1.0 + m_e) / (1.0 - m_e)) :
        (m_kind == Kind::Hyperbolic)
        ? SqRt ((m_e + 1.0) / (m_e - 1.0))
        : Inf<double>;

      // Mean Motion and its analogues: Surprisingly, defined in all cased:
      m_n =
        (m_kind == Kind::Circular     || m_kind == Kind::Elliptic   ||
         m_kind == Kind::RectElliptic || m_kind == Kind::Hyperbolic ||
         m_kind == Kind::RectHyperbolic)
        ? SqRt (Mu / Cube(m_a)) * 1.0_rad     // NB: Need "Abs" for *Hyperb*
        :
        (m_kind == Kind::Parabolic)
        ? 2.0_rad * SqRt (Mu / Cube(m_p))    
        : 6.0_rad * SqRt (Mu / Cube(1.0_km)); // RectParabolic

      // But the Orbital Period is only available for periodic orbits, other-
      // wise it is +oo:
      m_P =
        (m_kind == Kind::Circular     || m_kind == Kind::Elliptic   ||
         m_kind == Kind::RectElliptic)
        ? TWO_PI / m_n
        : Time(Inf<double>);

      //---------------------------------------------------------------------//
      // The position in orbit @ "a_t0": "m_E0", "m_f0":                     //
      //---------------------------------------------------------------------//
      // We initially characterise it by the Eccentric Anomaly (or its variants)
      // rather than by True Anomaly, because the latter is not informative for
      // RectLin orbits,   and is less accurate for highly-eccentrical Elliptic
      // orbits. However, it is undefined for Circular orbits!
      //
      // Also need the following:    Are we approaching the PeriFocus or moving
      // away from it?
      bool movingAway = IsPos(DotProd(a_pos, a_vel));

      switch (m_kind)
      {
      //------------------------//
      case Kind::Elliptic:
      case Kind::RectElliptic:
      //-----------------------//
      {
        // Eccentric Anomaly @ "m_t0":
        assert(0 < m_e && m_e <= 1);
        m_E0 = Angle(ACos((1.0 - double(r / m_a)) / m_e));
        // m_E0 in [0 .. Pi]

        if (!movingAway)
          m_E0 = TWO_PI - m_E0;
        assert(0.0_rad <= m_E0 && m_E0 < TWO_PI);

        // True Anomaly @ "a_t0":
        if (m_kind == Kind::Elliptic)
        {
          m_f0 =
            (LIKELY(m_E0 != PI)
            ? Angle(2.0 * ATan(m_S * Tan(m_E0 / 2.0)))
            : PI);
          // Initially in (-Pi, Pi]; then:
          m_f0 = To2Pi(m_f0);
        }
        else
          m_f0 = PI;  // RectElliptic
        assert(0.0_rad <  m_f0 && m_f0 < TWO_PI);

        // Mean Anomaly @ "a_t0" (Kepler's Equation):
        m_M0 = m_E0 - Angle(m_e * Sin(m_E0));
        assert(0.0_rad <= m_E0 && m_E0 < TWO_PI);
        break;
      }
      //-----------------------//
      case Kind::Hyperbolic:
      case Kind::RectHyperbolic:
      //-----------------------//
      {
        assert(m_e > 1.0);
        // We must have chE0 >= 1; otherwise, correct a rounding error:
        double chE0 = std::max((1.0 + double(r / m_a)) / m_e, 1.0);

        // XXX: In this case, it would be better to keep "m_E0" dim-less rather
        // than treat it as "Angle_rad", but we have to do the latter  for type
        // uniformity:
        m_E0 = Angle(ACosH(chE0));
        if (!movingAway)
          m_E0 = - m_E0;

        // True Anomaly @ "a_t0": This is really an Angle:
        if (m_kind == Kind::Hyperbolic)
        {
          m_f0 = Angle(2.0 * ATan(m_S * TanH(double(m_E0 / 2.0))));
          assert(-PI < m_f0 && m_f0 < PI);
        }
        else
          // RectHyperbolic:
          m_f0 = PI;

        // The Mean Anomaly / Kepler's Equation analogue (works for the Rect-
        // Hyperbolic orbit as well):
        m_M0 = Angle(m_e * SinH(double(m_E0))) - m_E0;
        break;
      }
      //-----------------------//
      case Kind::Parabolic:
      //-----------------------//
      {
        assert(m_e == 1.0);
        // m_E0 = tan(m_f0/2) is the argument of Barker's Equation. Again, this
        // not really an Angle_rad,  but we have to treat it as such due to our
        // type system constraints:
        m_E0 = Angle(SqRt(double(r / m_q) - 1.0));
        if (!movingAway)
          m_E0 = - m_E0;

        // True Anomaly @ "a_t0": This is really an Angle:
        m_f0 = Angle(2.0 * ATan(double(m_E0)));
        assert(-PI < m_f0 && m_f0 < PI);

        // The Mean Anomaly analogue  (consistent with "m_n" defined above):
        m_M0 = m_E0 * (1.0 + Sqr(double(m_E0)) / 3.0);
        break;
      }
      //-----------------------//
      case Kind::RectParabolic:
      //-----------------------//
      {
        // Interestingly, this case is different from "Parabolic":
        m_E0 = Angle(SqRt(2.0 * double(r / 1.0_km)));
        if (!movingAway)
          m_E0 = - m_E0;
        m_M0 = Angle(Cube(double(m_E0)));  // (consistent with "m_n" above)
        m_f0 = PI;
        break;
      }
      //-----------------------//
      case Kind::Circular:
      //-----------------------//
      {
        assert(!isRect);       // Of course...
        // In this case, m_M0 = m_E0 = m_f0;
        // set them all to the Longitude in Orbit "m_l0" (from the Ascneding
        // Node, or if I=0 or I=Pi, from the X axis):
        m_M0 = m_E0 = m_f0 = m_l0;
        break;
      }
      default:
        assert(false);
      }
      //---------------------------------------------------------------------//
      // Argument of the PeriFocus:                                          //
      //---------------------------------------------------------------------//
      if (!(isRect || m_kind == Kind::Circular))
      {
        // Generic Case: The PeriFocus exists:
        assert(m_kind == Kind::Elliptic  || m_kind == Kind::Parabolic ||
               m_kind == Kind::Hyperbolic);

        if (!(Iis0 || IisPi))
        {
          // The Ascneding Node exists:
          assert(IsFinite(m_l0) && IsFinite(m_f0));
          m_omega = To2Pi(m_l0  - m_f0);
        }
        else
        {
          // The Orbit is in the XY plane, no Ascending Node: "m_omega" is set
          // to the Longitude of the PeriFocus from X of "COS":
          assert(IsZero(m_Omega));
          // And "m_l0" is the longitude in Orbit from X, so:
          m_omega = To2Pi((Iis0 ? m_l0 : -m_l0) - m_f0);
        }
      }
      else
        // For RectiLinear or Circular Orbits, "m_omega" is totally undefined,
        // and is set to 0:
        m_omega = 0.0_rad;

      //---------------------------------------------------------------------//
      // The PeriFocus Passage Time: the last one before "a_t0":             //
      //---------------------------------------------------------------------//
      // In all cases, we have "m_M0" and "m_n" defined:
      assert(IsPos(m_n));
      Time tau = m_M0 / m_n;
      m_t0     = a_t0;
      m_T      = a_t0 - tau;

      //---------------------------------------------------------------------//
      // Construct the Rotation Matrices: "COS" <-> Orbital CoOrdSystem:     //
      //---------------------------------------------------------------------//
      MkTMs();

#     ifndef NDEBUG
      //---------------------------------------------------------------------//
      // End-to-End Check:                                                   //
      //---------------------------------------------------------------------//
      // Re-Construct the "a_pos" and "a_vel":
      Vector3D<LenK, COS, B> pos0;
      Vector3D<VelK, COS, B> vel0;
      GetPV(m_f0,  &pos0, &vel0);
      assert(pos0.ApproxEquals(a_pos));
      assert(vel0.ApproxEquals(a_vel));
#     endif
    }

    //=======================================================================//
    // "MkEllipticOrbit":                                                    //
    //=======================================================================//
    static TwoBodyOrbit<COS,B> MkEllipticOrbit
    (
      LenK    a_a,
      double  a_e,
      Angle   a_I,
      Angle   a_Omega,
      Angle   a_omega,
      TDB     a_T
    )
    {
      assert(0.0     <  a_e && a_e <  1.0 && IsPos(a_a) &&
             0.0_rad <= a_I && a_I <= PI);

      TwoBodyOrbit<COS,B> eo;

      // Main Orbital Elements:
      eo.m_kind  = Kind::Elliptic;
      eo.m_a     = a_a;
      eo.m_e     = a_e;
      // I, Omega, omega and l0 are deferred...
      eo.m_T     = a_T;

      // Derived Orbital Elements:
      eo.m_q     = eo.m_a * (1.0 - eo.m_e);
      eo.m_Q     = eo.m_a * (1.0 + eo.m_e);
      eo.m_p     = eo.m_a * (1.0 - Sqr(eo.m_e));
      eo.m_h     = SqRt(Mu   * eo.m_p);
      eo.m_S     = SqRt((1.0 + eo.m_e)   /  (1.0 - eo.m_e));
      eo.m_n     = SqRt(Mu / Cube(eo.m_a)) * 1.0_rad;
      eo.m_P     = TWO_PI  / eo.m_n;
      // "m_t0" is set to the PeriFocus Time, and then:
      eo.m_t0    = eo.m_T;
      eo.m_M0    = 0.0_rad;
      eo.m_f0    = 0.0_rad;

      // Orbit Orientation:
      bool Iis0  = a_I < Angle (1e-6);
      bool IisPi = a_I > (1.0 - 1e-6) * PI;

      if (UNLIKELY(Iis0))
      {
        // Ascending Node is undefined, then special conventions apply:
        eo.m_I     = 0.0_rad;
        eo.m_Omega = 0.0_rad;
        // "a_omega" is then interpreted as the PeriFocus Longitude from the
        // X axis of the "COS":
        eo.m_l0    = To2Pi(a_omega);
        // Argument of the PeriFocus is counted from the X axis:
        eo.m_omega = To2Pi(eo.m_l0 - eo.m_f0);  // "To2Pi" not really needed
      }
      else
      if (UNLIKELY(IisPi))
      {
        // Similar to the case I=0 above:
        eo.m_I     = PI;
        eo.m_Omega = 0.0_rad;
        eo.m_l0    = To2Pi(a_omega);
        eo.m_omega = To2Pi(- eo.m_l0 - eo.m_f0);
      }
      else
      {
        // Generic Case:
        eo.m_I     = a_I;
        eo.m_Omega = To2Pi(a_Omega);
        eo.m_omega = To2Pi(a_omega);
        eo.m_l0    = eo.m_omega;
      }

      // Finally, construct the Rotation Matrices:
      eo.MkTMs();
      return eo;
    }

    //=======================================================================//
    // "TimeSincePeriFocus":                                                 //
    //=======================================================================//
    // NB: Because the arg is the True Anomaly "a_f", this function is UNDEFINED
    // for RectiLinear orbits (where f==Pi==const):
    //
    constexpr Time TimeSincePeriFocus(Angle a_f)
    {
      assert(m_kind != Kind::RectElliptic && m_kind != Kind::RectParabolic &&
             m_kind != Kind::RectHyperbolic);

      // Mean Motion is always defined:
      assert(IsPos(m_n));

      switch (m_kind)
      {
      case Kind::Circular:
        // The PeriFocus is notional; M=E=f=l:
        return a_f / m_n;

      case Kind::Elliptic:
      case Kind::RectElliptic:
      {
        // Get the Eccentric Anomaly:
        assert(m_S > 0.0);
        Angle E =
          (a_f != PI)
          ? Angle(2.0 * ATan(Tan(a_f / 2.0) / m_S))
          : PI;
        E = To2Pi(E);

        // Mean Anomaly (Kepler's Equation):
        Angle  M = To2Pi(E - Angle(m_e * Sin(E)));
        return M / m_n;
      }

      case Kind::Hyperbolic:
      {
        assert(m_e > 1.0 && m_S > 0.0);

        // Eccentric Anomaly analogue:
        // Unlike the Ctor, here we do not need "E" and "M" to be Angles:
        double E = 2.0 * ATanH(Tan(a_f / 2.0) / m_S);

        // Mean Anomaly (Kepler's Equation analogue):
        double M = m_e * SinH(E) - E;
        return (M * 1.0_rad) / m_n;
      }

      case Kind::Parabolic:
      {
        assert(m_e == 1.0);

        // Eccentric Anomaly analogue:
        double E = Tan(a_f / 2.0);

        // Mean Anomaly analogue (Barker's Equation):
        double M = E * (Sqr(E) / 3.0);
        return (M * 1.0_rad)   / m_n;
      }

      default:
        assert(false);
      }
      __builtin_unreachable();
    }

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    constexpr LenK         a()        const { return m_a;  }
    constexpr double       e()        const { return m_e;  }
    constexpr Angle        I()        const { return m_I;  }
    constexpr Angle        Omega()    const { return m_Omega; }
    constexpr Angle        omega()    const { return m_omega; }
    constexpr TDB          T()        const { return m_T;  }
    constexpr LenK         q()        const { return m_q;  }
    constexpr LenK         Q()        const { return m_Q;  }
    constexpr LenK         p()        const { return m_p;  }
    constexpr AngVel       n()        const { return m_n;  }
    constexpr Time         P()        const { return m_P;  }
    constexpr Angle        M0()       const { return m_M0; }
    constexpr Angle        l0()       const { return m_l0; }
    constexpr TDB          t0()       const { return m_t0; }

    constexpr Mtx33 const& GetTM()    const { return m_TM;    }
    constexpr Mtx33 const& GetInvTM() const { return m_invTM; }

    //=======================================================================//
    // Utils:                                                                //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Radius-Vector as a function of True Anomaly:                          //
    //-----------------------------------------------------------------------//
    constexpr LenK GetR(Angle a_f) const
      { return m_p / (1.0 + m_e * Cos(a_f)); }

    //-----------------------------------------------------------------------//
    // Position and Velocity Vectors as functions of True Anomaly:           //
    //-----------------------------------------------------------------------//
    constexpr void GetPV
    (
      Angle a_f,
      Vector3D<LenK, COS, B>* a_pos, // Non-NULL
      Vector3D<VelK, COS, B>* a_vel  // May be NULL
    )
    const
    {
      double cosf = Cos(a_f);
      double sinf = Sin(a_f);

      //---------------------------------------------------------------------//
      // Position:                                                           //
      //---------------------------------------------------------------------//
      assert(a_pos != nullptr);

      // In-Orbit Co-Ords:
      LenK r =    GetR(a_f);
      LenK        inOrbP[3] { r * cosf, r * sinf, 0.0_km };
      // Perform the rotation:
      m_TM.MVMult(inOrbP, a_pos->GetArr());

      //---------------------------------------------------------------------//
      // Velocity:                                                           //
      //---------------------------------------------------------------------//
      if (a_vel != nullptr)
      {
        auto r2  = Sqr(r);

        // Total Velocity, from the Energy Integral:
        DEBUG_ONLY(VelK V = SqRt(Mu * (2.0 / GetR(a_f) - 1.0 / m_a));)

        // Angular Velocity, from the 2nd Kepler's Law:
        // r^2 * fDot = const = h:
        auto   fDot   = m_h / r2;

        // In-Orbit Velocity Components: rDot and (r * fDot):
        // "rDot"   is @ ( cos(f), sin(f));
        // "VTrans" is @ (-sin(f), cos(f)):
        VelK   rDot   = r2 / m_p * m_e * sinf * fDot;
        VelK   VTrans = r  * fDot;
        assert((Sqr(rDot)  + Sqr(VTrans)).ApproxEquals(Sqr(V)));

        // In-Orbit Velocity Vector:
        VelK inOrbV[3]
        {
          rDot * cosf - VTrans * sinf,
          rDot * sinf + VTrans * cosf,
          VelK(0.0)
        };
        // Perform the rotation:
        m_TM.MVMult(inOrbV, a_vel->GetArr());
      }
    }

    //-----------------------------------------------------------------------//
    // Angle Normalisation to the Range [0 .. 2*Pi):                         //
    //-----------------------------------------------------------------------//
    constexpr static Angle To2Pi(Angle a_phi)
    {
      // Normally, "a_phi" would not be far from the target range, so adjust it
      // in the following simple way:
      while (IsNeg(a_phi))
        a_phi += TWO_PI;
      while (a_phi >= TWO_PI)
        a_phi -= TWO_PI;
      assert(0.0_rad <= a_phi && a_phi < TWO_PI);
      return a_phi;
    }

  private:
    //=======================================================================//
    // "MkTMs":                                                              //
    //=======================================================================//
    void MkTMs()
    {
      // Let [X1, Y1, Z1] be an intermediate Co-Ords system with:
      // X1: PeriFocus -> Ascending Node of the Orbital Plane
      // Z1: Same as Z'   (normal to the Orbital Plane, CCW Motion from Z1 end)
      // Y1: [X1, Y1, Z1] is a right-oriented frame.
      // Then r = TM1 * r1, where:
      Mtx33 TM1;

      double cosOm = Cos(m_Omega);
      double sinOm = Sin(m_Omega);
      double cosI  = Cos(m_I);
      double sinI  = Sin(m_I);

      // XXX: The following sign factor "s" accounts for the "RetroGrade" motion
      // when I > Pi/2:
      double s = (m_I <= PI_2) ? 1.0 : -1.0;

      TM1(0,0) = cosOm; TM1(0,1) = -s*cosI*sinOm; TM1(0,2) =  s*sinI*sinOm;
      TM1(1,0) = sinOm; TM1(1,1) =  s*cosI*cosOm; TM1(1,2) = -s*sinI*cosOm;
      TM1(2,0) = 0.0;   TM1(2,1) =    sinI;       TM1(2,2) =    cosI;

      // Then r1 =       R3(-omega) * r', and thus
      //      r  = TM1 * R3(-omega) * r':
      m_TM    = TM1 * Mtx33::MkR3(-m_omega);
      m_invTM = m_TM.Transpose();
    }
  };
}
// End namespace SpaceBallistics
