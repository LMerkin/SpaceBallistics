// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/ME/TrConeSpherSegm.hpp":               //
//                  Geometric Objects as "Mechanical Elements"               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/TrConeSpherSegm.h"
#include "SpaceBallistics/ME/MechElement.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TrCone" Non-Default Ctor:                                              //
  //=========================================================================//
  // See "RotationShell" for the meaning of "{cos|sin}{A|P}":
  //
  template<LVSC    LVSCKind>
  constexpr TrCone<LVSCKind>::TrCone
  (
    Len     a_xu,        // Upper Base Center (Upper=Larger-X axis end)
    Len     a_yu,        //
    Len     a_zu,        //
    double  a_cosA,      // cos(alpha), >  0
    double  a_sinA,      // sin(alpha), >= 0
    double  a_cosP,      // cos(psi),   any sign
    double  a_sinP,      // sin(psi),   any sign
    Len     a_du,        // Upper (Larger-X)  Base Diameter
    Len     a_dl,        // Lower (Smaller-X) Base Diameter
    Len     a_h,         // Height
    Density a_rho,       // Propellant Density (0 if no Propellant)
    Mass    a_empty_mass // May be UnKnown
  )
  {
    //-----------------------------------------------------------------------//
    // Over-All Geometry:                                                    //
    //-----------------------------------------------------------------------//
    assert( IsPos(a_h)    && !IsNeg(a_du)  && !IsNeg(a_dl) &&
          !(IsZero(a_du)  && IsZero(a_dl)) && a_cosA > 0.0 &&
            a_sinA >= 0.0 && !IsNeg(a_rho) && !IsNeg(a_empty_mass));
    // XXX: Further geometry checks are performed by "RS::Init"...

    // The over-all sizes:
    Len   R     = a_du / 2.0;   // Upper base radius
    Len   r     = a_dl / 2.0;   // Lower base radius
    m_deltaR    = R - r;        // Any sign
    Len2  h2    = Sqr(a_h);

    // Memoised coeffs of the Cardano formula (used in "PropLevelOfVol"):
    m_clD       = h2  / Pi<double>;
    m_clVol     = 3.0 * m_deltaR * m_clD;
    m_rh        = r * a_h;
    m_rh3       = Cube(m_rh);
    m_PiR2      = Pi<double> * Sqr(R);

    //-----------------------------------------------------------------------//
    // Parent Classes Initialisation:                                        //
    //-----------------------------------------------------------------------//
    // For Optimisation:
    Len    s    = SqRt(Sqr(m_deltaR) + h2);
    double a    = double(m_deltaR /  a_h);
    double a2   = Sqr(a);
    double a3   = a2 * a;
    double a4   = Sqr(a2);
    Len2   R2   = Sqr(R);
    Len2   r2   = Sqr(r);
    Len3   r3   = r2 * r;
    Len4   r4   = Sqr(r2);

    // Side Surface Area and Nominal Enclosed Volume:
    Area   sideSurfArea = Pi<double>     * s   * (R  + r);
    Vol    enclVol      = Pi<double>/3.0 * a_h * (R2 + R * r + r2);
    m_baseArea          = enclVol   /a_h;  // Meaningful for a Cylinder only

    // "Intrinsic" "empty" MoIs:
    Len4   JE0  = Pi<double> * h2 * s * (R / 2.0 + r / 6.0);
    Len4   JE1  = Pi_4<double>    * s * (R + r)  * (R2  + r2);
    Len3   KE   = Pi<double>/3.0  * s *  a_h * (2.0 * R + r );

    // Coeffs of "intrtinsic" MoIs with Propellant:
    m_JP05  = Pi<double> / 5.0 * a2;
    m_JP04  = Pi_2<double>     * a  * r;
    m_JP03  = Pi<double> / 3.0 * r2;

    m_JP15  = Pi<double> /20.0 * a4;
    m_JP14  = Pi_4<double>     * a3 * r;
    m_JP13  = Pi_2<double>     * a2 * r2;
    m_JP12  = Pi_2<double>     * a  * r3;
    m_JP11  = Pi_4<double>     * r4;

    m_KP4   = Pi_4<double>     * a2;
    m_KP3   = TwoPi<double> / 3.0   * a * r;
    m_KP2   = Pi_2<double>     * r2;

    // Similar Coeffs for JP0, JP1 and KP Rates ("Dots"):
    m_JP0D4 = 5.0 * m_JP05;
    m_JP0D3 = 4.0 * m_JP04;
    m_JP0D2 = 3.0 * m_JP03;

    m_JP1D4 = 5.0 * m_JP15;
    m_JP1D3 = 4.0 * m_JP14;
    m_JP1D2 = 3.0 * m_JP13;
    m_JP1D1 = 2.0 * m_JP12;
    m_JP1D0 =       m_JP11;

    m_KPD3  = 4.0 * m_KP4;
    m_KPD2  = 3.0 * m_KP3;
    m_KPD1  = 2.0 * m_KP2;

    // Initialise the Parent Classes' Flds: NB: For (xu,yu,zu), IsUp=true:
    RS::Init
    (
      sideSurfArea,  enclVol, a_empty_mass,  a_cosA, a_sinA, a_cosP, a_sinP,
      a_xu, a_yu, a_zu, true, a_h, JE0, JE1, KE,     a_rho
    );
  }

  //=========================================================================//
  // "TrCone" Non-Default Ctors, Simple Cases (TrCone/Cylinder w/ OX axis):  //
  //=========================================================================//
  // "TrCone" with the rotation axis coinsiding with Ox: yu=zu=alpha=0, psi
  // is irrelevant and is also set to 0:
  //
  template<LVSC    LVSCKind>
  constexpr TrCone<LVSCKind>::TrCone
  (
    Len     a_xu,         // Upper (Larger-X)  Base Center
    Len     a_du,         // Upper (Larger-X)  Base Diameter
    Len     a_dl,         // Lower (Smaller-X) Base Diameter
    Len     a_h,          // Height
    Density a_rho,        // 0 may be OK
    Mass    a_empty_mass  // May be UnKnown
  )
  : TrCone(a_xu, 0.0_m, 0.0_m, 1.0, 0.0, 1.0, 0.0, a_du, a_dl, a_h, a_rho,
           a_empty_mass)
  {}

  // As above, but with dU=dL, ie a Cylinder:
  template<LVSC    LVSCKind>
  constexpr TrCone<LVSCKind>::TrCone
  (
    Len     a_xu,         // Upper (Larger-X)  Base Center
    Len     a_d,          // Diameter of both  Bases
    Len     a_h,          // Height
    Density a_rho,        // 0 may be OK
    Mass    a_empty_mass  // May be UnKnown
  )
  : TrCone(a_xu, 0.0_m, 0.0_m, 1.0, 0.0, 1.0, 0.0, a_d, a_d, a_h, a_rho,
           a_empty_mass)
  {}

  //=========================================================================//
  // "TrCone" Propellant Volume -> Propellant Level:                         //
  //=========================================================================//
  // Returns (Level, LevelDot):
  //
  template<LVSC LVSCKind>
  constexpr std::pair<Len,Vel> TrCone<LVSCKind>::PropLevelOfVol
    (Vol a_v, VolRate a_v_dot)
  const
  {
    assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= ME::GetEnclVol());
    bool isFull = (a_v  == ME::GetEnclVol());

    Len l   (NaN<double>);
    Vel lDot(NaN<double>);

    if (IsZero(m_deltaR))
    {
      // R==r, ie a Cylinder:
      assert(IsPos(m_baseArea));
      l    = isFull  ? RS::GetHeight() : (a_v / m_baseArea);
      lDot = a_v_dot / m_baseArea;
    }
    else
    {
      // A (possibly Truncated) Cone indeed:
      if (isFull)
      {
        l    = RS::GetHeight();
        lDot =
          // XXX: Beware of R==0:
          LIKELY(IsPos(m_PiR2))
          ? a_v_dot / m_PiR2
          :
          IsPos(a_v_dot)
          ? Vel( Inf<double>)
          :
          IsNeg(a_v_dot)
          ? Vel(-Inf<double>)
          : Vel(0.0);
      }
      else
      {
        // General case: Solving a cubic equation by the Cardano formula;
        // since Vol'(l) > 0 globally, there is only 1 real root:
        // m_deltaR != 0:
        auto x = CbRt(m_clVol * a_v  + m_rh3);
        assert(!IsNeg(x));

        l      = (x - m_rh) / m_deltaR;
        lDot   =
          // It is possible that x==0, but only if a_v==0 && r==0:
          LIKELY(IsPos(x))
          ? m_clD / Sqr(x) * a_v_dot
          :
          IsPos(a_v_dot)
          ? Vel (Inf<double>)
          :
          IsNeg(a_v_dot)
          ? Vel(-Inf<double>)
          : Vel(0.0);
      }
      assert(!(IsNeg(l) || IsPos(lDot)));
      return std::make_pair(l, lDot);
    }
    assert(IsFinite(l) && IsFinite(lDot) && !(IsNeg(l) || IsPos(lDot)));
    return std::make_pair(l, lDot);
  }

  //=========================================================================//
  // "TrCone" MoI Components for the Propellant of Given Level:              //
  //=========================================================================//
  // These functions are exactly the same for "TrCone" and "SpherSegm", yet we
  // cannot factor them out without creating an intermediate parent class.  So
  // use a macro to generate them:
  //
# ifdef  MK_PROP_MOI_COMPS
# undef  MK_PROP_MOI_COMPS
# endif
# define MK_PROP_MOI_COMPS(ClassName)  \
  template<LVSC LVSCKind>     \
  constexpr void ClassName<LVSCKind>::PropMoIComps \
  ( \
    Len       a_l,       \
    Vel       a_l_dot,   \
    Len5*     a_jp0,     \
    Len5*     a_jp1,     \
    Len4*     a_kp,      \
    Len5Rate* a_jp0_dot, \
    Len5Rate* a_jp1_dot, \
    Len4Rate* a_kp_dot   \
  ) \
  const \
  { \
    assert(!(IsNeg(a_l)         || IsPos(a_l_dot))      && \
           a_l       <= RS::GetHeight()                 && \
           a_jp0     != nullptr && a_jp1     != nullptr && \
           a_kp      != nullptr && a_jp0_dot != nullptr && \
           a_jp1_dot != nullptr && a_kp_dot  != nullptr);  \
    \
    Len2 l2 = Sqr (a_l); \
    Len3 l3 = l2 * a_l;  \
    Len4 l4 = Sqr (l2);  \
    Len5 l5 = l4 * a_l;  \
    \
    /* MoI Components: */   \
    *a_jp0 =  m_JP05 * l5 + m_JP04 * l4 + m_JP03 * l3; \
    *a_jp1 =  m_JP15 * l5 + m_JP14 * l4 + m_JP13 * l3 + m_JP12 * l2 + \
              m_JP11 * a_l; \
    *a_kp  =  m_KP4  * l4 + m_KP3  * l3 + m_KP2  * l2; \
    \
    /* MoI Component's "Dots": */ \
    *a_jp0_dot = (m_JP0D4 * l4  + m_JP0D3  * l3 + m_JP0D2 * l2)  * a_l_dot; \
    *a_jp1_dot = (m_JP1D4 * l4  + m_JP1D3  * l3 + m_JP1D2 * l2   +          \
                  m_JP1D1 * a_l + m_JP1D0) * a_l_dot;                       \
    *a_kp_dot  = (m_KPD3  * l3  + m_KPD2   * l2 + m_KPD1  * a_l) * a_l_dot; \
    \
    assert(!(IsNeg(*a_jp0)     || IsNeg(*a_jp1)     || IsNeg(*a_kp)    ||   \
             IsPos(*a_jp0_dot) || IsPos(*a_jp1_dot) || IsPos(*a_kp_dot)));  \
  }

  // Apply the above macro:
  MK_PROP_MOI_COMPS(TrCone)

  //=========================================================================//
  // "SpherSegm" Non-Default Ctor:                                           //
  //=========================================================================//
  template<LVSC LVSCKind>
  constexpr SpherSegm<LVSCKind>::SpherSegm
  (
    bool    a_facing_up,
    Len     a_xb,          // Base center
    Len     a_yb,          //
    Len     a_zb,          //
    double  a_cosA,        // As in "RotationShell"
    double  a_sinA,        //
    double  a_cosP,        //
    double  a_sinP,        //
    Len     a_d,           // Base diameter
    Len     a_h,           // Height
    Density a_rho,         // Propellant Density (may be 0)
    Mass    a_empty_mass   // May be UnKnown
  )
  {
    //-----------------------------------------------------------------------//
    // Over-All Geometry:                                                    //
    //-----------------------------------------------------------------------//
    assert(IsPos(a_d)   &&  IsPos(a_h) && a_cosA > 0.0 && a_sinA >= 0.0 &&
          !IsNeg(a_rho) && !IsNeg(a_empty_mass));
    // XXX: Furher geometry checks are performed by "RS::Init"...

    Len  r        = a_d / 2.0;                    // Base   radius
    assert(a_h   <= r * TolFact);                 // Important!
    a_h           = std::min(a_h, r);
    m_R           = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius
    m_R3          = Cube(m_R);
    m_cR3         = (Pi<double> / 3.0) * m_R3;
    m_facingUp    = a_facing_up;

    //---------------------------------------------------------------------//
    // Parent Classes Initialisation:                                      //
    //---------------------------------------------------------------------//
    // Side Surface Area and Notional Volume:
    Len2 h2     = Sqr(a_h);
    Len3 h3     = a_h * h2;
    Len2 R2     = Sqr(m_R);

    // Side Surface Area and Nominal Enclosed Volume:
    Area sideSurfArea = TwoPi<double> * m_R *  a_h;
    Vol  enclVol      =    Pi<double> * h2  * (m_R - a_h / 3.0);

    // "Intrinsic" "empty" MoIs (NB: they do not depend on "FacingUp",
    // because JE0 is wrt the rotation axis "xi"; furthermore, the Up- and
    // Low-facing segments  are mutually-symmetric wrt any axis orthogonal
    // to "xi", hence "JE1" must also be the same; for "KE",  it is just a
    // lucky coincidence due to f(xi)*sqrt(1+f'(xi)^2)=R=const):
    //
    Len4 JE0 = TwoPi<double> / 3.0 * m_R * h3;
    Len4 JE1 = m_R * enclVol;
    Len3 KE  = Pi<double> * m_R * h2;

    // Coeffs of "intrtinsic" MoIs with Propellant. Unlike the "empty" ones
    // above, these coeffs do depend on the Segment orientation:
    Len Rmh  = m_R - a_h;
    assert(!IsNeg(Rmh));
    Len TRmh = m_R + Rmh;

    m_JP05   = -Pi<double>   / 5.0;
    m_JP04   =  Pi_2<double> * (a_facing_up ? -Rmh : m_R);
    m_JP03   =  a_facing_up
                ? Pi<double> / 3.0 * TRmh * a_h
                : Len2(0.0);

    m_JP15   =  Pi<double>   / 20.0;
    m_JP14   =  Pi_4<double> * (a_facing_up ? Rmh : -m_R);
    m_JP13   =  Pi<double> *
                (a_facing_up
                 ? R2 / 3.0 - m_R * a_h + h2 / 2.0
                 : R2 / 3.0);
    m_JP12   =  a_facing_up
                ? -Pi_2<double> * Rmh * TRmh * a_h
                : Len3(0.0);
    m_JP11   =  a_facing_up
                ? Pi_4<double>  * Sqr(TRmh)  * h2
                : Len4(0.0);

    m_KP4    = -Pi_4<double>;
    m_KP3    =  (TwoPi<double> / 3.0) * (a_facing_up  ? -Rmh : m_R);
    m_KP2    =  a_facing_up
                ? Pi_2<double>  * TRmh * a_h
                : Len2(0.0);

    // Similar Coeffs for JP0, JP1 and KP Rates ("Dots"):
    m_JP0D4  = 5.0 * m_JP05;
    m_JP0D3  = 4.0 * m_JP04;
    m_JP0D2  = 3.0 * m_JP03;

    m_JP1D4  = 5.0 * m_JP15;
    m_JP1D3  = 4.0 * m_JP14;
    m_JP1D2  = 3.0 * m_JP13;
    m_JP1D1  = 2.0 * m_JP12;
    m_JP1D0  =       m_JP11;

    m_KPD3   = 4.0 * m_KP4;
    m_KPD2   = 3.0 * m_KP3;
    m_KPD1   = 2.0 * m_KP2;

    // Initialise the Parent Classes' Flds:
    // NB: (xb,yb,zb) is the Base Center, so for it, BaseIsUp = !IsFacingUp:
    RS::Init
    (
      sideSurfArea, enclVol, a_empty_mass, a_cosA,  a_sinA, a_cosP, a_sinP,
      a_xb, a_yb, a_zb, !a_facing_up, a_h,    JE0, JE1, KE, a_rho
    );
  }

  //=========================================================================//
  // "SperSegm" Non-Default Ctors, Simple Cases:                             //
  //=========================================================================//
  // Rotation axis coinciding with OX, so alpha=0; in this case, we can assume
  // psi=0 as well:
  //
  template<LVSC LVSCKind>
  constexpr SpherSegm<LVSCKind>::SpherSegm
  (
    bool    a_facing_up,
    Len     a_xb,         // Base center X co-ord
    Len     a_d,          // Base diameter
    Len     a_h,          // Height
    Density a_rho,        // 0 may be OK
    Mass    a_empty_mass  // May be UnKnown
  )
  : SpherSegm(a_facing_up, a_xb,  0.0_m, 0.0_m, 1.0, 0.0, 1.0, 0.0,
              a_d,    a_h, a_rho, a_empty_mass)
  {}

  // As above, but with d/2 = h, ie a HemiSpehere:
  //
  template<LVSC LVSCKind>
  constexpr SpherSegm<LVSCKind>::SpherSegm
  (
    bool    a_facing_up,
    Len     a_xb,         // Base center X co-ord
    Len     a_d,          // Base diameter
    Density a_rho,        // 0 may be OK
    Mass    a_empty_mass  // May be UnKnown
  )
  : SpherSegm(a_facing_up,    a_xb,  0.0_m, 0.0_m, 1.0, 0.0, 1.0, 0.0,
              a_d, a_d / 2.0, a_rho, a_empty_mass)
  {}

  //=========================================================================//
  // "SpherSegm" Propellant Volume -> Propellant Level:                      //
  //=========================================================================//
  template<LVSC LVSCKind>
  constexpr std::pair<Len,Vel> SpherSegm<LVSCKind>::PropLevelOfVol
    (Vol a_v, VolRate a_v_dot)
  const
  {
    assert(!(IsNeg(a_v) || IsPos(a_v_dot)));
    return
      m_facingUp
      ? PropLevelOfVolUp (a_v, a_v_dot)
      : PropLevelOfVolLow(a_v, a_v_dot);
  }

  //=========================================================================//
  // "PropLevelOfVol" for the Low-Facing "SpherSegm":                        //
  //=========================================================================//
  template<LVSC LVSCKind>
  constexpr std::pair<Len,Vel> SpherSegm<LVSCKind>::PropLevelOfVolLow
    (Vol a_v, VolRate a_v_dot)
  const
  {
    assert(!IsNeg(a_v) && a_v <= ME::GetEnclVol());

    // However, "a_v_dot" may be of any sign, as this function may be called
    // either from "PropLevelOfVol" or from "PropLevelOfVolUp"...
    //
    Len l   (NaN<double>);
    Vel lDot(NaN<double>);

    if (UNLIKELY(IsZero(a_v)))
    {
      // NB: This case corresponds to infinite "lDot", be careful of its sign.
      // When a_v_dot==0, "lDot" is obviously 0:
      l    = Len(0.0);
      lDot = IsPos(a_v_dot)
             ? Vel( Inf<double>) :
             IsNeg(a_v_dot)
             ? Vel(-Inf<double>)
             : Vel(0.0);
    }
    else
    {
      // Generic Case: a_v > 0:
      // Solving the equation
      //   x^2*(3-x) = tv,
      // where
      //   "x"  is the Propellant level relative to "R": x=l/R, 0 <= x  <= 1;
      //   "tv" is the Volume/(Pi/3*R^3),                       0 <= tv <= 2;
      // by using Halley's Method (to avoid dealing with complex roots in the
      // Cardano formula):
      // x=0.5 is a reasonble initial approximation (XXX: We may use "h"  for
      // a more accurate initial point, but this is not required yet):
      //
      double x    = 0.5;
      double tv   = double(a_v / m_cR3);
      assert(0.0 <= tv && tv < 2.0 * TolFact);
      tv = std::min(2.0,  tv);           // Enforce the boundary

      // For safety, restrict the number of iterations:
      constexpr int N = 100;
      int           i = 0;
      for (; i < N; ++i)
      {
        double x2 = Sqr(x);
        double x3 = x2 * x;
        double x4 = Sqr(x2);
        double dx =
          x * (x - 2.0)   * (x3 - 3.0 * x2 + tv) /
          (2.0 * x4 - 8.0 *  x3 + 9.0 * x2 + tv * (1.0 - x));

        // Iterative update:
        x -= dx;

        // Exit condition:
        if (UNLIKELY(Abs(dx) < DefaultTol<double>))
          break;
      }
      // If we got here w/o achieving the required precision, it's an error:
      assert(i < N);

      // If all is fine: To prevent rounding errors, enforce the boundaries:
      x = std::max(std::min(x, 1.0), 0.0);

      // So the Level:
      l = x * m_R;

      // And finally, "lDot":
      // Derive it from the above cubic equation:
      // 3*x*(2-x) * x_dot = tv_dot ;   the LHS always exists:
      // BEWARE of x==0 and other boundary cases:
      if (LIKELY(x > 0.0))
      {
        assert(IsPos(m_cR3));
        auto tvDot = a_v_dot / m_cR3;
        auto xDot  = tvDot   / (3.0 * x * (2.0 - x));
        lDot       = xDot    * m_R;
      }
      else
        lDot       =
          IsPos(a_v_dot)
          ? Vel( Inf<double>) :
          IsNeg(a_v_dot)
          ? Vel(-Inf<double>)
          : Vel(0.0);
    }
    assert(IsFinite(l) && !IsNeg(l)); // But "lDot" can be any...
    return std::make_pair(l, lDot);
  }

  //=========================================================================//
  // "PropLevelOfVol" for the Up-Facing "SpherSegm":                         //
  //=========================================================================//
  // Using the invariant
  // V_up(l) + V_low(h-l) = EnclVol:
  //
  template<LVSC LVSCKind>
  constexpr std::pair<Len,Vel> SpherSegm<LVSCKind>::PropLevelOfVolUp
    (Vol a_v, VolRate a_v_dot)
  const
  {
    assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= ME::GetEnclVol());
    auto lowRes = PropLevelOfVolLow(ME::GetEnclVol() - a_v, - a_v_dot);
    Len l       = RS::GetHeight() - lowRes.first;
    Vel lDot    =                 - lowRes.second;
    assert(!(IsNeg(l) || IsPos(lDot)));
    return std::make_pair(l, lDot);
  }

  //=========================================================================//
  // "SpherSegm" MoI Components for the Propellant of Given Level:           //
  //=========================================================================//
  // Exacyly same function as for "TrCone", so apply the above macro:
  MK_PROP_MOI_COMPS(SpherSegm)
# undef MK_PROP_MOI_COMPS
}
// End namespace SpaceBallistics
