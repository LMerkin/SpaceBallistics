// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/EllipticOrbit.hpp":                 //
//                  Elliptic Keplerian (Two-Body) Motion                     //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/PhysForces/BodyData.hpp"
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
    Angle  m_l0;     // Longitude from the Ascending Node @ m_t0
    Angle  m_M0;     // Mean Anomaly @ m_t0

    // Default Ctor is meaningless:
    EllipticOrbit() = delete;

  public:
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
      auto h      = angMom.EuclNorm();
      assert(IsPos(h));   // FIXME: Rectilinear ellipses are not allowed

      //  "k" is the unit vector orthogonal to the Orbital Plane:
      Vector3D<double, COS, B> k = angMom / h;

      // Orbit Inclination:
      m_I     = Angle(ACos(k.z()));
      assert(0.0_rad <= m_I && m_I <= Angle(Pi<double>));

      // FIXME: For the moment, we disallow the degenerate case I==0 or I==Pi,
      // since Omega and omega are not well-defined then:
      assert(0.0_rad <  m_I && m_I <  Angle(Pi<double>));

      // ProGrade or RetroGrade motion?
      bool isProGrade = m_I <= Angle(Pi_2<double>);

      // Ascending Node Longitude:
      // k.x  = [+-]sin(I) * sin(Omega), k.y = [-+]sin(I) * cos(Omega),
      // sin(I) >= 0 always; FIXME: Is the case I==Pi/2 treated correctly
      // by default? Also, I=0 and I=Pi need to be verified...
      m_Omega =
        Angle(isProGrade
              ? ATan2( k.x(), -k.y())
              : ATan2(-k.x(),  k.y()));
      if (IsNeg(m_Omega))
        m_Omega += Angle(TwoPi<double>);

      // Semi-Major Axis, from the Energy Integral:
      // V^2  = Mu * (2/r - 1/a) :
      m_a     = 1.0 / (2.0 / r - Sqr(V) / Mu);
      assert(IsPos(m_a));

      // Focal Param ("p") and the Eccentricity:
      m_p     = Sqr(h) / Mu;
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
      m_P     = Angle(TwoPi<double>)  / m_n;

      // "f": True Anomaly @ "a_t0":
      // p/r = 1 + e * cos(f);
      // "f" is defined by cos(f) up to the sign; the latter is determined
      // by the (pos,vel) angle:
      Angle f = Angle(ACos((double(m_p / r) - 1.0) / m_e));
      if (IsNeg(DotProd(a_pos, a_vel)))
        f = -f;

      // Longitude in orbit from the Ascending Node:
      double cosl0 =
        double((a_pos.x() * Cos(double(m_Omega))  +
                a_pos.y() * Sin(double(m_Omega))) / r);
      m_l0 = Angle(ACos(cosl0));      // In [0 .. Pi] initially

      // Possibly adjust "u" in the generic case (FIXME: under the assumption
      // that we are in the generic case, ie I != 0 and I != Pi):
      if (IsNeg(a_pos.z()))
        // Then, obviously, Pi < u < 2*Pi, so we must update it but preserve
        // the "cosU":
        m_l0 = Angle(TwoPi<double>) - m_l0;
      assert(0.0_rad <= m_l0 && m_l0 < Angle(TwoPi<double>));

      // We can now get the Argument of the PeriFocus:
      m_omega = m_l0 - f;
      if (IsNeg(m_omega))
        m_omega += Angle(TwoPi<double>);
      else
      if (m_omega >= Angle(TwoPi<double>))
        m_omega -= Angle(TwoPi<double>);
      assert(0.0_rad <= m_omega && m_omega < Angle(TwoPi<double>));

      m_pi  = m_Omega + m_omega;
      if (m_pi >= Angle(TwoPi<double>))
        m_pi -=   Angle(TwoPi<double>);
      assert(0.0_rad <= m_pi && m_pi < Angle(TwoPi<double>));

      // Finally, the PeriFocus Passage Time -- the last one before "a_t0":
      // FIXME: it is undefined if e=0, but that is not allowed yet:
      // Eccentric Anomaly:
      Angle E = Angle(2.0 * ATan(Tan(double(f) / 2.0) / m_S));
      assert  (-Angle(Pi<double>) < E && E < Angle(Pi<double>));

      // Mean Anomaly (Kepler's Equation):
      m_M0 = E - Angle(m_e * Sin(double(E)));
      if (IsNeg(m_M0))
        m_M0 += Angle(TwoPi<double>);
      assert(0.0_rad <= m_M0 && m_M0 < Angle(TwoPi<double>));

      // Time from PeriFocus and the PeriFocus Passage Time:
      Time tau = m_M0 / m_n;
      m_t0     = a_t0;
      m_T      = a_t0 - tau;
    }

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    constexpr LenK   a()     const { return m_a; }
    constexpr double e()     const { return m_e; }
    constexpr Angle  I()     const { return m_I; }
    constexpr Angle  Omega() const { return m_Omega; }
    constexpr Angle  omega() const { return m_omega; }
    constexpr TDB    T()     const { return m_T; }
    constexpr LenK   q()     const { return m_q; }
    constexpr LenK   Q()     const { return m_Q; }
    constexpr LenK   p()     const { return m_p; }
    constexpr AngVel n()     const { return m_n; }
    constexpr Angle  M0()    const { return m_M0; }
    constexpr Angle  l0()    const { return m_l0; }
    constexpr TDB    t0()    const { return m_t0; }
  };
}
// End namespace SpaceBallistics
