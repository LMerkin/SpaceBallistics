// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/ME/ToricSegms.hpp":                  //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/MechElement.hpp"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "ToricSegm" Class:                                                      //
  //=========================================================================//
  // Describes the Upper or Lower part of a Toric Shell; "r" and "h"  are the
  // Radius and Height of the Segment (similar to those of SpherSegm, but un-
  // like the "SpherSegm", the "ToricSegm" is obtained by rotation of a circu-
  // lar arc around an axis which is parallel to the arc axis, but not coinci-
  // ding with the latter). "R" is the Minor Radius of the Torus (derived from
  // "r" and "h" similar to "SpherSegm"), and "Q" is  the Major Radius of the
  // Torus, that is, the distance between the aforemention arc axis and the ro-
  // tation axis; we may assume r <= R <= Q.  Up-/Low-Facing orientation is si-
  // milar to that of a "SpherSegm":
  //
  class ToricSegm final: public RotationBody<ToricSegm>
  {
  private:
    using RotationBody<ToricSegm>::Tol;
    using RotationBody<ToricSegm>::TolFact;

    //=======================================================================//
    // Data Flds: Toric Segment's Geometry:                                  //
    //=======================================================================//
    Len     m_Q;        // Major Radius of the Toric Segm (NOT Torus!)
    Len     m_R;        // Minor Radius of the Torus      (not just Segm!)
    bool    m_facingUp; // Orientation
    double  m_q2;       // Memoised Coeffs
    double  m_q2L;      //
    Len     m_twoH;     //
    Len2    m_h2;       //
    Vol     m_NV;       //
    Len4    m_pQR3;     //
    Len5    m_pQR4;     //
    Len5    m_JPL0;     // MoI components for Full-Propellant in Low-Facing Sgm
    Len5    m_JPL1;     //
    Len4    m_KPL;      //

    // Default Ctor is deleted:
    ToricSegm() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr ToricSegm
    (
      bool        a_facing_up,
      Len         a_xb,          // Over-All   Segment Base Center
      Len         a_yb,          //
      Len         a_zb,          //
      double      a_alpha,       // Dimension-less, in (-Pi/2 .. Pi/2)
      Len         a_d,           // Cross-Section Base Diameter
      Len         a_h,           // Cross-Section Height
      Len         a_D,           // Over-All   Segment Diameter
      Density     a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass        a_empty_mass = UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // Geometry:                                                           //
      //---------------------------------------------------------------------//
      assert(IsPos(a_d)   && IsPos(a_h) && Abs(a_alpha) < Pi_2<double> &&
            !IsNeg(a_rho) && IsPos(a_D) && !IsNeg(a_empty_mass));

      Len r       =  a_d / 2.0;                   // Cross-Section Base Radius
      m_Q         =  a_D / 2.0 - r;               // Major Segment      Radius
      assert(a_h  <= r * TolFact && m_Q > r);
      a_h         = std::min (a_h,  r);
      m_R         = (Sqr(r) / a_h + a_h) / 2.0;   // Torus Minor Radius
      m_facingUp  = a_facing_up;
      a_h         = std::min (a_h,  m_R);         // For safety

      // We also assume that m_Q >= R, although this is not strictly necessary:
      assert(m_Q >= m_R && m_R >= a_h);

      // Memoised Coeffs:
      Len2 R2     = Sqr(m_R);
      m_q2        = Sqr(double(m_Q / m_R));
      assert(m_q2 >= 1.0);
      m_q2L       = m_q2 + 0.75;
      m_twoH      = 2.0 * a_h;
      m_h2        = Sqr  (a_h);
      m_NV        = TwoPi<double>    * Sqr(m_R) * m_Q;
      Len3 pQR2   = Pi<double> * m_Q * R2;
      m_pQR3      = pQR2   * m_R;
      m_pQR4      = m_pQR3 * m_R;

      // Pre-Compute the MoI Components for Full Propellant Load for Low-Facing
      // Segm, and the Nominal Enclosed Volume. The "Dots" are ignored:
      Len5Rate JPL0Dot, JPL1Dot;
      Len4Rate KPLDot;
      Vol      enclVol;
      VolRate  enclVolRate;
      PropMoICompsLow
      (
        a_h,      Vel(0.0),
        &m_JPL0,  &m_JPL1,  &m_KPL,
        &JPL0Dot, &JPL1Dot, &KPLDot,
        &enclVol, &enclVolRate
      );
      assert
        (IsPos (m_JPL0)  && IsPos (m_JPL1)  && IsPos (m_KPL)   &&
         IsPos (enclVol) && IsZero(JPL0Dot) && IsZero(JPL1Dot) &&
         IsZero(KPLDot)  && IsZero(enclVolRate));

      //---------------------------------------------------------------------//
      // Parent Classes Initialisation:                                      //
      //---------------------------------------------------------------------//
      // Side Surface Area:
      double x   = double(a_h / m_R);
      assert(0.0 < x && x <= 1.0);
      double cx  = 1.0 - x;
      double acx = std::acos(cx);  // XXX: Again, not a "constexpr" in CLang...
      double s   = SqRt(x * (2.0 - x));
      Area   sideSurfArea = (4.0 * Pi<double> * acx) * m_R * m_Q;
      assert(IsPos(sideSurfArea));

      // "Intrinsic" "empty" MoI components:
      Len4 JE0 = (2.0 * m_pQR3) *
                 (m_facingUp
                  ? (3.0 - 4.0*x + 2.0*Sqr(x)) * acx - 3.0  * cx * s
                  :  3.0                       * acx - (x + 3.0) * s);
      Len4 JE1 = m_pQR3 *  ((2.0 * m_q2 + 3.0) * acx - 3.0  * cx * s);
      Len3 KE  = (4.0 * pQR2) * (m_facingUp  ? s - cx * acx : acx - s);

      // Initialise the Parent Classes' Flds:
      // NB: (xb,yb,zb) is the Base Center, so BaseIsUp = !IsFacingUp:
      RotationBody<ToricSegm>::Init
      (
        sideSurfArea,  enclVol,    a_empty_mass,
        a_alpha, a_xb, a_yb, a_zb, !a_facing_up, a_h,
        JE0,     JE1,  KE,   a_rho
      );
    }

    //=======================================================================//
    // Non-Default Ctor, Simple Cases:
    //=======================================================================//
    // Rotation axis coinciding with OX:
    //
    constexpr ToricSegm
    (
      bool        a_facing_up,
      Len         a_xb,       // Over-All Segment Base Center (X-CoOrd)
      Len         a_d,        // Cross-Section Base Diameter
      Len         a_h,        // Cross-Section Height
      Len         a_D,        // Over-All Segment Diameter
      Density     a_rho,      // 0 may be OK (if holds no Propellant)
      Mass        a_empty_mass = UnKnownMass
    )
    : ToricSegm(a_facing_up, a_xb, 0.0_m, 0.0_m, 0.0, a_d, a_h, a_D,
                a_rho,       a_empty_mass)
    {}

    // As above, but with d/2 = h, ie the Cross-Section is a HemiSpehere:
    //
    constexpr ToricSegm
    (
      bool        a_facing_up,
      Len         a_xb,       // Over-All Segment Base Center (X-CoOrd)
      Len         a_d,        // Cross-Section Base Diameter
      Len         a_D,        // Over-All Segment Diameter
      Density     a_rho,      // 0 may be OK (if holds no Propellant)
      Mass        a_empty_mass = UnKnownMass
    )
    : ToricSegm(a_facing_up, a_xb, 0.0_m, 0.0_m, 0.0, a_d, a_d / 2.0, a_D,
                a_rho,       a_empty_mass)
    {}

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const
    {
      assert(!(IsNeg(a_v) || IsPos(a_v_dot)));
      return
        m_facingUp
        ? PropLevelOfVolUp (a_v, a_v_dot)
        : PropLevelOfVolLow(a_v, a_v_dot);
    }

  private:
    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolLow
      (Vol a_v, VolRate a_v_dot) const
    {
      assert(!IsNeg(a_v));
      // But "a_v_dot" may be of any sign in different use cases...

      // Solving the equation
      //   arccos(z) - z * SqRt(1-z^2) = y
      // where
      //   z = 1-x,
      //   "x" is the Propellant level relative to "R": x=l/R, 0 <= x <= 1;
      //   "y" is the Volume/(NV = 2*Pi*R^2*Q),                0 <= y <= Pi/2;
      // by using Halley's Method;
      // x=z=1/2 is a reasonble initial approximation:
      //
      double z  = 0.5;
      double y  = double(a_v / m_NV);
      assert(0.0 <= y && y < Pi_2<double> * TolFact);
      y = std::min (y, Pi_2<double>);       // Enforce the upper boundary

      // For safety, restrict the number of iterations:
      constexpr int N = 100;
      int           i = 0;
      for (; i < N; ++i)
      {
        assert(-Tol < z && z < TolFact);
        z = std::min(std::max(z, 0.0), 1.0);

        double z2   = Sqr(z);
        double cz2  = 1.0 - z2;
        assert(cz2 >= 0.0);
        double s    = SqRt(cz2);
        double dy   = y - std::acos(z);
        double dz   =
          2.0 * cz2 * (dy + z * s) / ((4.0 - 3.0 * z2) * s + z * dy);

        // Iterative update:
        z -= dz;

        // Exit condition:
        if (UNLIKELY(Abs(dz) < Tol))
          break;
      }
      // If we got here w/o achieving the required precision, it's an error:
      assert(i < N);

      // If all is fine: To prevent rounding errors, enforce the boundaries:
      z = std::min(std::max(z, 0.0), 1.0);

      // So the Propellant level:
      Len l = (1.0 - z) * m_R;
      assert(!IsNeg(l));

      // And "lDot" from the equation solved:
      // -2 * sqrt(1-z^2) * z_dot = y_dot :
      //
      auto yDot = a_v_dot    / m_NV;
      Vel  lDot = m_R * yDot / (2.0 * SqRt(1.0 - Sqr(z)));  // Infinite if z=1
      assert(!IsPos(lDot));

      return std::make_pair(l, lDot);
    }

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    // Using the invariant
    // V_up(l) + V_low(h-l) = EnclVol:
    // Same implementation as for "SpherSegm::PropLevelOfVolUp":
    //
    constexpr std::pair<Len,Vel> PropLevelOfVolUp
      (Vol a_v, VolRate a_v_dot) const
    {
      assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= GetEnclVol());

      auto lowRes = PropLevelOfVolLow(GetEnclVol() - a_v, - a_v_dot);
      Len l       = GetHeight() - lowRes.first;
      Vel lDot    = -             lowRes.second;
      assert(!IsNeg(l) || IsPos(lDot));
      return std::make_pair (l, lDot);
    }

  public:
    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr void PropMoIComps
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const
    {
      assert(!(IsNeg(a_l) || IsPos(a_l_dot)) && a_l <= GetHeight());
      if (m_facingUp)
        PropMoICompsUp
          (a_l, a_l_dot, a_jp0, a_jp1, a_kp, a_jp0_dot, a_jp1_dot, a_kp_dot);
      else
        PropMoICompsLow
          (a_l, a_l_dot, a_jp0, a_jp1, a_kp, a_jp0_dot, a_jp1_dot, a_kp_dot,
           nullptr,      nullptr);
    }

    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    // XXX: This method is also called from the "ToricSegm" Ctor, when the par-
    // ent classes are not yet initialised, so it MUST NOT call any methods of
    // "RotationBody":
    //
    constexpr void PropMoICompsLow
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot,
      Vol*      a_vp,       // May be NULL
      VolRate*  a_vp_dot    // Ditto
    )
    const
    {
      // NB: "a_l_dot" may be for any sign in different use cases of this func:
      assert(!IsNeg(a_l)                                  &&
             a_jp0     != nullptr && a_jp1     != nullptr &&
             a_kp      != nullptr && a_jp0_dot != nullptr &&
             a_jp1_dot != nullptr && a_kp_dot  != nullptr);

      // In this case, the formulas are relatively simple, so use the direct
      // computations:
      double x   = double(a_l / m_R);
      assert(-Tol < x &&  x <= TolFact);
      x = std::min(std::max(x, 0.0), 1.0);

      double cx  = 1.0 - x;
      double acx = std::acos(cx);  // XXX: Again, not a "constexpr" in Clang...
      double x2  = Sqr(x);
      double x3  = x2 * x;
      double s   = SqRt(x * (2.0 - x));

      // MoI Components:
      *a_jp0     =
        m_pQR4   * (2.5   * acx + s * (x3 - x2/3.0 - (5.0/6.0) * x - 2.5));
      *a_jp1     =
        m_pQR4   * (m_q2L * acx - s * cx * (m_q2L  + x - 0.5 * x2));
      *a_kp      =
         m_pQR3  * (2.0   * acx - s * (x + 1.0) * (2.0 - (4.0/3.0) * x));

      // And their "Dots":
      auto xDot  = a_l_dot / m_R;
      auto jkDot = 4.0 * x * s * xDot;
      *a_jp0_dot = m_pQR4  * x * jkDot;
      *a_jp1_dot = m_pQR4  * (2.0 * s * (m_q2L + 2.0 * x - x2 - 0.75)) * xDot;
      *a_kp_dot  = m_pQR3  * jkDot;

      if (a_vp != nullptr)
        *a_vp   = m_NV * (acx - cx  * s);

      if (a_vp_dot != nullptr)
        *a_vp_dot   = m_NV  *  (2.0 * s) * xDot;
    }

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    // Here there is no optional rerturn values "VP" and "VPDot":
    //
    constexpr void PropMoICompsUp
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const
    {
      assert(!(IsNeg(a_l)         || IsPos(a_l_dot))      &&
             a_jp0     != nullptr && a_jp1     != nullptr &&
             a_kp      != nullptr && a_jp0_dot != nullptr &&
             a_jp1_dot != nullptr && a_kp_dot  != nullptr);

      // In this case, the formulas are more involved,  so we first compute the
      // MoI components of the corresp complementary Segment, take a difference
      // wrt a Fully-Propellant-Loaded Segment, and apply a shift:
      //
      // MoI components for a complementary segment  can be obtained using the
      // same formulas as for a Low-Facing one (for which the Propellant is ad-
      // jacent to the Pole, not to the Bottom Plane):
      Len5     JC0,    JC1;
      Len4     KC;
      Vol      VC;
      Len5Rate JC0Dot, JC1Dot;
      Len4Rate KCDot;
      VolRate  VCDot;

      PropMoICompsLow
      (
        GetHeight() - a_l, -a_l_dot,
        &JC0,    &JC1,     &KC,
        &JC0Dot, &JC1Dot,  &KCDot,
        &VC,     &VCDot
      );
      // Subtract the above vals from the MoI components for a Fully-Loaded
      // Low-Facing Segment:
      Len5 JP0 = m_JPL0 - JC0;
      Len5 JP1 = m_JPL1 - JC1;
      Len4 KP  = m_KPL  - KC;
      Vol  VP  = GetEnclVol() - VC;
      assert(IsPos(JP0) && IsPos(JP1) && IsPos(KP) && IsPos(VP));

      // "JP0" is now almost what we need:
      // for our  required segm,  Base    @ 0,     Surface @ l,
      // for what we got so far,  Surface @ (h-l), Base    @ h;
      // the orientation does not matter, but the distance to the OEta axis is;
      // so move it in the Low direction of Xi by "h":
      *a_jp0     =   JP0    + m_h2 * VP    - m_twoH * KP;
      *a_jp0_dot = -(JC0Dot + m_h2 * VCDot - m_twoH * KCDot);
      assert(IsPos(*a_jp0) && !IsPos(*a_jp0_dot));

      // "JP1" is taken as is, because the Eta co-ords are unaffected by Xi
      // shifts:
      *a_jp1     =  JP1;
      *a_jp1_dot = -JC1Dot;
      assert(IsPos(*a_jp1) && !IsPos(*a_jp1_dot));

      // "KP" is subject to a similar shift and sign inversion (due to mirror
      // symmetry wrt OEta axis):
      *a_kp      =   GetHeight() * VP    - KP;
      *a_kp_dot  = -(GetHeight() * VCDot - KCDot);
      assert(IsPos(*a_kp)  && !IsPos(*a_kp_dot));

      // All Done!
    }
  };

  //=========================================================================//
  // "DoubleCylinder" Class:                                                 //
  //=========================================================================//
  // Provides a "cylindrical torus" (a body obtained by rotating a rectangle
  // around an outside axis parallel to the cylinder's main axis):
  //
  class DoubleCylinder: public RotationBody<DoubleCylinder>
  {
  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Memoised Coeffs for Propellant MoI Components:
    Area    m_S;       // Base Area
    Len2    m_cJP0;
    Len2    m_dJP0;
    Len4    m_cJP1;
    Len2    m_cKP;
    Len2    m_dKP;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The General Case:                                                     //
    //-----------------------------------------------------------------------//
    constexpr DoubleCylinder
    (
      Len         a_xu,          // Upper  Base Center
      Len         a_yu,          //
      Len         a_zu,          //
      double      a_alpha,       // Dimension-less, in (-Pi/2 .. Pi/2)
      Len         a_D,           // Outer Diameter
      Len         a_d,           // Inner Diameter
      Len         a_h,           // Height
      Density     a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass        a_empty_mass = UnKnownMass
    )
    {
      assert(IsPos(a_D)  && IsPos(a_d)   && a_D > a_d     && IsPos(a_h) &&
             Abs(a_alpha) < Pi_2<double> && !IsNeg(a_rho) &&
             !IsNeg(a_empty_mass));

      // Geometry:
      Len   R  = a_D / 2.0;
      Len   r  = a_d / 2.0;
      Len2  R2 = Sqr(R);
      Len2  r2 = Sqr(r);
      Len2  h2 = Sqr(a_h);
      Len3  R3 = R * R2;
      Len3  r3 = r * r2;
      Len4  R4 = Sqr(R2);
      Len4  r4 = Sqr(r2);

      // Side Surface Area and Enclosed Volume:
      Area   sideSurfArea = Pi<double> * (a_D + a_d) * a_h;
      Vol    enclVol      = Pi<double> * (R2  - r2 ) * a_h;

      // "Intrinsic" "empty" MoI components:
      Len4 JE0 = sideSurfArea * h2  / 3.0;
      Len4 JE1 = Pi<double>   * a_h * (R3 + r3);
      Len3 KE  = Pi<double>   * h2  * (R  + r );

      // Memoised coeffs for the Propellant MoIs:
      m_S      = Pi<double>   * (R2 - r2);
      m_cJP0   = m_S / 3.0;
      m_dJP0   = 3.0 * m_cJP0;
      m_cJP1   = Pi_4<double> * (R4 - r4);
      m_cKP    = m_S / 2.0;
      m_dKP    = 2.0 * m_cKP;

      // Initialise the Parent Classes' Flds:
      // NB: (xu,yu,zu) is the Upper Base Center, so BaseIsUp = true here:
      RotationBody<DoubleCylinder>::Init
      (
        sideSurfArea,  enclVol,    a_empty_mass,
        a_alpha, a_xu, a_yu, a_zu, true,   a_h,
        JE0,     JE1,  KE,   a_rho
      );
    }

    //-----------------------------------------------------------------------//
    // Simple Case: Rotation Axis coincides with OX:                         //
    //-----------------------------------------------------------------------//
    constexpr DoubleCylinder
    (
      Len         a_xu,          // Upper  Base Center
      Len         a_D,           // Outer Diameter
      Len         a_d,           // Inner Diameter
      Len         a_h,           // Height
      Density     a_rho,         // Propellant Density (may be 0 if no Propelt)
      Mass        a_empty_mass = UnKnownMass
    )
    : DoubleCylinder
        (a_xu, 0.0_m, 0.0_m, 0.0, a_D, a_d, a_h, a_rho, a_empty_mass)
    {}

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const
    {
      assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= GetEnclVol());
      Len l    = a_v     /  m_S;
      Vel lDot = a_v_dot /  m_S;
      return std::make_pair(l, lDot);
    }

    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr void PropMoIComps
    (
      Len       a_l,
      Vel       a_l_dot,
      Len5*     a_jp0,
      Len5*     a_jp1,
      Len4*     a_kp,
      Len5Rate* a_jp0_dot,
      Len5Rate* a_jp1_dot,
      Len4Rate* a_kp_dot
    )
    const
    {
      assert(!(IsNeg(a_l) || IsPos(a_l_dot))   && a_l <= GetHeight() &&
             a_jp0     != nullptr && a_jp1     != nullptr  &&
             a_kp      != nullptr && a_jp0_dot != nullptr  &&
             a_jp1_dot != nullptr && a_kp_dot  != nullptr);

      Len2 l2 = Sqr(a_l);
      Len3 l3 = a_l * l2;

      *a_jp0     = m_cJP0 * l3;
      *a_jp1     = m_cJP1 * a_l;
      *a_kp      = m_cKP  * l2;

      *a_jp0_dot = m_dJP0 * l2  * a_l_dot;
      *a_jp1_dot = m_cJP1       * a_l_dot;
      *a_kp_dot  = m_dKP  * a_l * a_l_dot;
    }
  };
}
// End namespace SpaceBallistics
