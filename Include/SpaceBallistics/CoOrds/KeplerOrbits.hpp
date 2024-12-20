// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/KeplerOrbits.hpp":                 //
//              Kinematics of the Keplerian (Two-Body) Motion                //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/PhysForces/BodyData.hpp"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "EllipticOrbit" Class:                                                  //
  //=========================================================================//
  template<typename COS, Body B = Body::UNDEFINED>
  class EllipticOrbit
  {
  private:
    //=======================================================================//
    // Elliptic Keplerian Elements:                                          //
    //=======================================================================//
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
    Angle  m_pi;     // PeriFocus Longitude
    LenK   m_q;      // PeriFocus Distance
    LenK   m_Q;      // ApoFocus  Distance
    LenK   m_p;      // Focal Parameter
    double m_S;      // SqRt((1+e)/(1-e))
    AngVel m_n;      // Mean Motion (rad/sec !!!)
    Time   m_P;      // Orbital Period  (sec !!!)
    TDB    m_t0;     // Time of these Elements
    Angle  m_M0;     // Mean Anomaly @ m_t0
    Angle  m_l0;     // Longitude from the Ascending Node @ m_t0
    Angle  m_f0;     // True Anomaly @ m_t0
    decltype(Sqr(1.0_km) / 1.0_sec)
           m_h;      // Angular Momentum Constant

    // "TM": Matrix for Transforming "In-Orbit" Co-Ords [X', Y', Z'] into the
    // COS, where:
    // X' axis: from the Focus to the PeriFocus;
    // Z' axis: from the Focus, normal to the orbital plane, so that the orbital
    //          motion is Counter-Clock-Wise when viewed from the end of Z';
    // Y' axis: such that [X',Y',Z'] is the right-oriented frame:
    //
    Mtx33  m_TM;     // Rotation Matrix: R = TM * R'
    Mtx33  m_invTM;  // Inverse (= Transposed)  "TM"

    // Default Ctor is meaningless:
    EllipticOrbit() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor: From the Orbital Elements:                          //
    //=======================================================================//
    constexpr EllipticOrbit
    (
      LenK    a_a,
      double  a_e,      // In   [0,1)
      Angle   a_I,      // In   [0..Pi]
      Angle   a_Omega,  // Normalised to [0..2*Pi)
      Angle   a_omega,  // ditto
      TDB     a_T       // Any
    )
    : // Main Orbital Elements:
      m_a     (a_a),
      m_e     (a_e),
      m_I     (a_I),
      m_Omega (To2Pi(a_Omega)),
      m_omega (To2Pi(a_omega)),
      m_T     (a_T),
      // Derived Orbital Elements:
      m_pi    (To2Pi(m_Omega + m_omega)),
      m_q     (m_a * (1.0 - m_e)),
      m_Q     (m_a * (1.0 + m_e)),
      m_p     (m_a * (1.0 - Sqr(m_e))),
      m_S     (SqRt((1.0 + m_e) / (1.0 - m_e))),
      m_n     (SqRt(Mu / Cube(m_a)) * 1.0_rad),
      m_P     (TWO_PI  / m_n),
      m_t0    (m_T),   // Set to the PeriFocus Time
      m_M0    (0.0_rad),
      m_l0    (m_omega),
      m_f0    (0.0_rad),
      m_h     (SqRt(Mu * m_p))
    {
      assert
        (IsPos(m_a) && 0.0 <= m_e && m_e < 1.0 && 0.0_rad <= m_I && m_I <= PI);
      // Then construct the Rotation Matrices:
      MkTMs();
    }

    //=======================================================================//
    // Non-Default Ctor: From Pos and Vel:                                   //
    //=======================================================================//
    constexpr EllipticOrbit
    (
      TDB                           a_t0,
      Vector3D<LenK, COS, B> const& a_pos,
      Vector3D<VelK, COS, B> const& a_vel
    )
    {
      // Both vectors must be non-0:
      LenK r = LenK(a_pos);
      VelK V = VelK(a_vel);
      assert(IsPos(r) && IsPos(V));

      // Construct the Angular Momentum Vector:
      auto angMom = CrossProd(a_pos, a_vel);
      m_h         = angMom.EuclNorm();
      assert(IsPos(m_h));   // FIXME: Rectilinear ellipses are not allowed yet

      //  "k" is the unit vector orthogonal to the Orbital Plane:
      Vector3D<double, COS, B> k = angMom / m_h;

      // Orbit Inclination:
      m_I     = Angle(ACos(k.z()));
      assert(0.0_rad <= m_I && m_I <= PI);

      // FIXME: For the moment, we disallow the degenerate case I==0 or I==PI,
      // since Omega and omega are not well-defined then:
      assert(0.0_rad <  m_I && m_I <  PI);

      // ProGrade or RetroGrade motion?
      bool isProGrade = m_I <= PI_2;

      // Ascending Node Longitude:
      // k.x  = [+-]sin(I) * sin(Omega), k.y = [-+]sin(I) * cos(Omega),
      // sin(I) >= 0 always; FIXME: Is the case I==PI/2 treated correctly
      // by default? Also, I=0 and I=PI need to be verified...
      m_Omega =
        Angle(isProGrade
              ? ATan2( k.x(), -k.y())
              : ATan2(-k.x(),  k.y()));
      m_Omega = To2Pi(m_Omega);

      // Semi-Major Axis, from the Energy Integral:
      // V^2  = Mu * (2/r - 1/a) :
      m_a     = 1.0 / (2.0 / r - Sqr(V) / Mu);
      assert(IsPos(m_a));

      // Focal Param ("p") and the Eccentricity:
      m_p     = Sqr(m_h) / Mu;
      assert(IsPos(m_a) && m_p <= m_a);

      m_e     = SqRt(1.0 - double(m_p / m_a));
      assert(0.0 <= m_e && m_e < 1.0);
      // FIXME: We currently disallow the circular motion case:
      assert(0.0 <  m_e && m_e < 1.0);

      // Derived Elements:
      m_q     = m_a * (1.0 - m_e);
      m_Q     = m_a * (1.0 + m_e);
      m_S     = SqRt ((1.0 + m_e) / (1.0 - m_e));
      m_n     = SqRt (Mu / Cube(m_a)) * 1.0_rad;
      m_P     = TWO_PI   / m_n;

      // "f": True Anomaly @ "a_t0":
      // p/r = 1 + e * cos(f);
      // "f" is defined by cos(f) up to the sign; the latter is determined
      // by the (pos,vel) angle:
      m_f0 = Angle(ACos((double(m_p / r) - 1.0) / m_e));
      if (IsNeg(DotProd(a_pos, a_vel)))
        m_f0 = -m_f0;
      m_f0   = To2Pi(m_f0);

      // Longitude in orbit from the Ascending Node:
      double cosl0 =
        double((a_pos.x() * Cos(double(m_Omega))  +
                a_pos.y() * Sin(double(m_Omega))) / r);
      m_l0 = Angle(ACos(cosl0));      // In [0 .. PI] initially

      // Possibly adjust "u" in the generic case (FIXME: under the assumption
      // that we are in the generic case, ie I != 0 and I != PI):
      if (IsNeg(a_pos.z()))
        // Then, obviously, PI < l0 < 2*PI, so we must update "l0" but preserve
        // the "cosl0":
        m_l0 = TWO_PI - m_l0;
      assert(0.0_rad <= m_l0 && m_l0 < TWO_PI);

      // We can now get the Argument of the PeriFocus:
      m_omega = To2Pi(m_l0    - m_f0);
      m_pi    = To2Pi(m_Omega + m_omega);

      // Finally, the PeriFocus Passage Time -- the last one before "a_t0":
      // FIXME: it is undefined if e=0, but that is not allowed yet:
      // "tau" is the Time since the PeriFocus as a function of "m_f0":
      auto [tau, M] = TimeSincePeriFocus(m_f0);
      m_t0 = a_t0;
      m_T  = a_t0 - tau;
      m_M0 = M;

      // And also construct the Rotation Matrices:
      MkTMs();

      // End-to-End Check:
      // Re-Construct the "a_pos" and "a_vel":
#     ifndef NDEBUG
      Vector3D<LenK, COS, B> pos0;
      Vector3D<VelK, COS, B> vel0;
      GetPV(m_f0,  &pos0, &vel0);
      assert(pos0.ApproxEquals(a_pos));
      assert(vel0.ApproxEquals(a_vel));
#     endif
    }

    //=======================================================================//
    // "TimeSincePeriFocus":                                                 //
    //=======================================================================//
    // Returns the Time and the curr MeanAnomaly:
    //
    constexpr std::pair<Time, Angle> TimeSincePeriFocus(Angle a_f)
    {
      // Just to make sure:
      a_f = To2Pi(a_f);

      // Eccentric Anomaly:
      Angle E =
        (a_f != PI)
        ? Angle(2.0 * ATan(Tan(double(a_f) / 2.0) / m_S))
        : PI;
      E = To2Pi(E);

      // Mean Anomaly (Kepler's Equation):
      Angle M = To2Pi(E - Angle(m_e * Sin(double(E))));

      // Time from PeriFocus and the PeriFocus Passage Time:
      return std::make_pair(M / m_n, M);
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
      { return m_p / (1.0 + m_e * Cos(double(a_f))); }

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
      double cosf   = Cos(double(a_f));
      double sinf   = Sin(double(a_f));

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

      double cosOm = Cos(double(m_Omega));
      double sinOm = Sin(double(m_Omega));
      double cosI  = Cos(double(m_I));
      double sinI  = Sin(double(m_I));

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
