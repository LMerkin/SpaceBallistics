// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/ME/TrConeSpherSegm.hpp":               //
//                  Geometric Objects as "Mechanical Elements"               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/MechElement.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TrCone" Class: Truncated Cone Shell (w/o and with Propellant):         //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only).
  // Geometry: Upper Base Diameter "du" (for Larger-X) and Lower Base Diameter
  //           "dl" (for Smaller-X), height "h"; may have  du<=>dl;
  //           either "du" or "dl" may be be 0 (ie full cone), but  not both.
  // Location: (xu,yu,zu) -- co-ords of the Upper (Larger-X) Base center; XXX:
  //           either "zu" or "yu" must be 0. The cone axis is lying in the XY
  //           or XZ plane, resp., at the angle  "alpha" to the positive direc-
  //           tion of the OX axis; we assume that |alpha| < Pi/2;  obviously,
  //           the "du" base diameter corresponds to the base center (xu,yu)
  //           (if zu=0, ie XY plane) or (xu,zu) (if yu=0, ie XZ plane); the
  //           "dl" diameter corresponds to the other (Smaller-X) end of the
  //           cone axis segment;
  //           when alpha=0, the axis of the Spherical Segment coincides with
  //           OX, which is the most common case.
  // Mass:     If the "mass" param is set, its value is used  for the element
  //           mass which is in this case final. If this param is 0, the mass
  //           is set under the assumtion of SutfDensity=1, and it is NOT fi-
  //           nal yet (needs to be set later via "ProRateMass").
  // Return value: The resulting "MechElement" object:
  //
  template<LVSC LVSCKind>
  class TrCone final: public RotationShell<LVSCKind, TrCone<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, TrCone<LVSCKind>>;

    //=======================================================================//
    // Data Flds: Truncated Cone's Geometry:                                 //
    //=======================================================================//
    // Memoise pre-computed coeffs of the Cardano formula in "LevelOfVol":
    Len       m_deltaR;
    Len2      m_clD;
    Len3      m_clVol;
    Len2      m_rh;
    Len6      m_rh3;

    // Coeffs for computation of "intrinsic" MoI params (J0, J1, K) as polynom-
    // ial functions (of degress 5 for JP0, JP1, and of degree 4 for KP) of the
    // Propellant Level, used in "PropMoIComps":
    double    m_JP05;     // coeff(JP0, l^5)
    Len       m_JP04;     // coeff(JP0, l^4)
    Len2      m_JP03;     // coeff(JP0, l^3)
    double    m_JP15;     // coeff(JP1, l^5)
    Len       m_JP14;     // coeff(JP1, l^4)
    Len2      m_JP13;     // coeff(JP1, l^3)
    Len3      m_JP12;     // coeff(JP1, l^2)
    Len4      m_JP11;     // coeff(JP1, l)
    double    m_KP4;      // coeff(KP,  l^4)
    Len       m_KP3;      // coeff(KP,  l^3)
    Len2      m_KP2;      // coeff(KP,  l^2)

    // And similar for JP0, JP1 and K Rates ("Dots"):
    double    m_JP0D4;
    Len       m_JP0D3;
    Len2      m_JP0D2;
    double    m_JP1D4;
    Len       m_JP1D3;
    Len2      m_JP1D2;
    Len3      m_JP1D1;
    Len4      m_JP1D0;
    double    m_KPD3;
    Len       m_KPD2;
    Len2      m_KPD1;

    // Default Ctor is deleted:
    TrCone() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr TrCone
    (
      Len       a_xu,       // Upper Base Center (Upper=Larger-X axis end)
      Len       a_yu,       //
      Len       a_zu,       //
      double    a_alpha,    // Dimension-less, in (-Pi/2 .. Pi/2)
      Len       a_du,       // Upper (Larger-X)  Base Diameter
      Len       a_dl,       // Lower (Smaller-X) Base Diameter
      Len       a_h,        // Height
      Density   a_rho,      // Propellant Density (0 if no Propellant)
      Mass      a_empty_mass = ME::UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // Over-All Geometry:                                                  //
      //---------------------------------------------------------------------//
      assert( IsPos(a_h)   && !IsNeg(a_du)  && !IsNeg(a_dl) &&
            !(IsZero(a_du) && IsZero(a_dl)) && Abs(a_alpha) < Pi_2<double> &&
            ! IsNeg(a_rho) && !IsNeg(a_empty_mass));

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

      //---------------------------------------------------------------------//
      // Parent Classes Initialisation:                                      //
      //---------------------------------------------------------------------//
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
        sideSurfArea,  enclVol, a_empty_mass,
        a_alpha, a_xu, a_yu,    a_zu,   true, a_h,
        JE0,     JE1,  KE,      a_rho
      );
    }

    //=======================================================================//
    // Non-Default Ctor, Simple Cases (TrCone/Cylinder with OX axis):        //
    //=======================================================================//
    // "TrCone" with the rotation axis coinsiding with Ox: yu=zu=alpha=0:
    //
    constexpr TrCone
    (
      Len         a_xu,       // Upper (Larger-X)  Base Center
      Len         a_du,       // Upper (Larger-X)  Base Diameter
      Len         a_dl,       // Lower (Smaller-X) Base Diameter
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = ME::UnKnownMass
    )
    : TrCone(a_xu, 0.0_m, 0.0_m, 0.0, a_du, a_dl, a_h, a_rho, a_empty_mass)
    {}

    // As above, but with dU=dL, ie a Cylinder:
    constexpr TrCone
    (
      Len         a_xu,       // Upper (Larger-X)  Base Center
      Len         a_d,        // Diameter of both  Bases
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = ME::UnKnownMass
    )
    : TrCone(a_xu, 0.0_m, 0.0_m, 0.0, a_d, a_d, a_h, a_rho, a_empty_mass)
    {}

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    // Returns (Level, LevelDot):
    constexpr std::pair<Len,Vel> PropLevelOfVol(Vol a_v, VolRate a_v_dot) const
    {
      assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= RS::GetEnclVol());

      if (IsZero(m_deltaR))
      {
        // R==r: The simplest and the most common case: A Cylinder:
        auto x    = RS::GetHeight() / RS::GetEnclVol();
        Len  l    = a_v     * x;
        Vel  lDot = a_v_dot * x;

        assert(!(IsNeg(l) || IsPos(lDot)));
        return std::make_pair(l, lDot);
      }
      else
      {
        // General case: Solving a cubic equation by the Cardano formula;
        // since Vol'(l) > 0 globally, there is only 1 real root:
        auto x    = CbRt(m_clVol * a_v  + m_rh3);
        Len  l    = (x - m_rh) / m_deltaR;
        Vel  lDot = m_clD      / Sqr(x) * a_v_dot;

        assert(!(IsNeg(l) || IsPos(lDot)));
        return std::make_pair(l, lDot);
      }
    }

    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    // This function is exacytly the same for "TrCone" and "SpherSegm", yet we
    // cannot factor it out without creating an untermediate parent class.  So
    // use a macro to generate it:
#   ifdef  MK_PROP_MOI_COMPS
#   undef  MK_PROP_MOI_COMPS
#   endif
#   define MK_PROP_MOI_COMPS  \
    constexpr void PropMoIComps \
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
    MK_PROP_MOI_COMPS
  };

  //=========================================================================//
  // "SpherSegm" Class:                                                      //
  //=========================================================================//
  // Similar to "TrCone", but for a Spherical Segment of the base diameter "d"
  // and the height (from the base plane to the pole) "h", where we assume
  // h <= d/2 (the equality corresponds to the case of a HemiSphere).
  // NB: this is a Shperical Segment,  NOT a Spherical Slice,  ie it always
  // contains a pole. The Spherical Slice would be a more general case and a
  // more close analogy of the Truncated Cone, but it would (arguably) have
  // almost no real applications.
  // Furthermore, "alpha" is the angle between the segment's axis and the po-
  // sitive direction of the X axis; we assume that |alpha| < Pi/2;
  // (xB, yB, zB) are the co-ords of the base center (NOT of the pole!), and
  // either "zB" or "yB" must be 0;  in the former case, the Segment's rotat-
  // ion axis is assumed to lie in the XY plane, in the latter -- in XZ.
  // The spherical segment may be facing "Up" (ie towards the positive direc-
  // tion of the OX axis), or otherise, as given by the "facing_up" flag.
  // When alpha=0, the axis of the Spherical Segment coincides with OX, which
  // is the most common case:
  //
  template<LVSC LVSCKind>
  class SpherSegm final: public RotationShell<LVSCKind, SpherSegm<LVSCKind>>
  {
  private:
    using ME = MechElement  <LVSCKind>;
    using RS = RotationShell<LVSCKind, SpherSegm<LVSCKind>>;

    //=======================================================================//
    // Data Flds: Spherical Segment's Geometry:                              //
    //=======================================================================//
    // Over-all sizes and Orientation:
    Len     m_R;          // Sphere  Radius
    Len3    m_R3;         // Sphere  Radius^3
    Len3    m_cR3;        // Pi*R^3 / 3 = HemiSphereVol / 2
    bool    m_facingUp;   // Segment Ortientation

    // Coeffs for computation of "intrinsic" MoI params (J0, J1, K) as polynom-
    // ial functions (of degress 5 for JP0, JP1, and of degree 4 for KP) of the
    // Propellant Level -- similar to "TrCone":
    //
    double  m_JP05;       // coeff(JP0, l^5)
    Len     m_JP04;       // coeff(JP0, l^4)
    Len2    m_JP03;       // coeff(JP0, l^3)
    double  m_JP15;       // coeff(JP1, l^5)
    Len     m_JP14;       // coeff(JP1, l^4)
    Len2    m_JP13;       // coeff(JP1, l^3)
    Len3    m_JP12;       // coeff(JP1, l^2)
    Len4    m_JP11;       // coeff(JP1, l)
    double  m_KP4;        // coeff(KP,  l^4)
    Len     m_KP3;        // coeff(KP,  l^3)
    Len2    m_KP2;        // coeff(KP,  l^2)

    // And similar for JP0, JP1 and K Rates ("Dots"):
    double  m_JP0D4;
    Len     m_JP0D3;
    Len2    m_JP0D2;
    double  m_JP1D4;
    Len     m_JP1D3;
    Len2    m_JP1D2;
    Len3    m_JP1D1;
    Len4    m_JP1D0;
    double  m_KPD3;
    Len     m_KPD2;
    Len2    m_KPD1;

    // Default Ctor is deleted:
    SpherSegm() = delete;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr SpherSegm
    (
      bool        a_facing_up,
      Len         a_xb,          // Base center
      Len         a_yb,          //
      Len         a_zb,          //
      double      a_alpha,       // Dimension-less, in (-Pi/2 .. Pi/2)
      Len         a_d,           // Base diameter
      Len         a_h,           // Height
      Density     a_rho,         // Propellant Density (may be 0)
      Mass        a_empty_mass = ME::UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // Over-All Geometry:                                                  //
      //---------------------------------------------------------------------//
      assert(IsPos(a_d)   &&  IsPos(a_h) && Abs(a_alpha) < Pi_2<double> &&
            !IsNeg(a_rho) && !IsNeg(a_empty_mass));

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
        sideSurfArea,   enclVol,    a_empty_mass,
        a_alpha,  a_xb, a_yb, a_zb, !a_facing_up, a_h,
        JE0,      JE1,  KE,   a_rho
      );
    }

    //=======================================================================//
    // Non-Default Ctor, Simple Cases:                                       //
    //=======================================================================//
    // Rotation axis coinciding with OX:
    //
    constexpr SpherSegm
    (
      bool        a_facing_up,
      Len         a_xb,       // Base centerr X co-ord
      Len         a_d,        // Base diameter
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = ME::UnKnownMass
    )
    : SpherSegm(a_facing_up, a_xb, 0.0_m, 0.0_m, 0.0, a_d, a_h, a_rho,
                a_empty_mass)
    {}

    // As above, but with d/2 = h, ie a HemiSpehere:
    //
    constexpr SpherSegm
    (
      bool        a_facing_up,
      Len         a_xb,       // Base center X co-ord
      Len         a_d,        // Base diameter
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = ME::UnKnownMass
    )
    : SpherSegm(a_facing_up, a_xb, 0.0_m, 0.0_m, 0.0, a_d, a_d / 2.0, a_rho,
                a_empty_mass)
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
    // For the Low-Facing "SpherSegm":                                       //
    //-----------------------------------------------------------------------//
    constexpr std::pair<Len,Vel> PropLevelOfVolLow
      (Vol a_v, VolRate a_v_dot) const
    {
      // Solving the equation
      //   x^2*(3-x) = tv,
      // where
      //   "x"  is the Propellant level relative to "R": x=l/R, 0 <= x  <= 1;
      //   "tv" is the Volume/(Pi/3*R^3),                       0 <= tv <= 2;
      // by using Halley's Method (to avoid dealing with complex roots in the
      // Cardano formula):
      // x=1/2 is a reasonble initial approximation (XXX: We may use "h"  for
      // a more accurate initial point, but this is not required yet):
      //
      assert(!IsNeg(a_v));
      // However, "a_v_dot" may be of any sign, as this function may be called
      // either from "PropLevelOfVol" or from "PropLevelOfVolUp"...

      double x   = 0.5;
      double tv  = double(a_v / m_cR3);
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
        if (UNLIKELY(Abs(dx) < Tol))
          break;
      }
      // If we got here w/o achieving the required precision, it's an error:
      assert(i < N);

      // If all is fine: To prevent rounding errors, enforce the boundaries:
      x = std::max(std::min(x, 1.0), 0.0);

      // The level:
      Len l = x * m_R;
      assert(!IsNeg(l));

      // And finally, "lDot":
      // Derive it from the above cubic equation:
      // 3*x*(2-x) * x_dot = tv_dot ;   the LHS always exists:
      //
      auto tvDot = a_v_dot / m_cR3;
      auto xDot  = tvDot   / (3.0 * x * (2.0 - x));
      Vel  lDot  = xDot    * m_R;

      return std::make_pair(l, lDot);
    }

    //-----------------------------------------------------------------------//
    // For the Up-Facing "SpherSegm":                                        //
    //-----------------------------------------------------------------------//
    // Using the invariant
    // V_up(l) + V_low(h-l) = EnclVol:
    //
    constexpr std::pair<Len,Vel> PropLevelOfVolUp
      (Vol a_v, VolRate a_v_dot) const
    {
      assert(!(IsNeg(a_v) || IsPos(a_v_dot)) && a_v <= RS::GetEnclVol());
      auto lowRes = PropLevelOfVolLow(RS::GetEnclVol() - a_v, - a_v_dot);
      Len l       = RS::GetHeight() - lowRes.first;
      Vel lDot    =                 - lowRes.second;
      assert(!(IsNeg(l) || IsPos(lDot)));
      return std::make_pair(l, lDot);
    }

  public:
    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    // Exacyly same function as for "TrCone", so apply the above macro:
    MK_PROP_MOI_COMPS
#   undef MK_PROP_MOI_COMPS
  };
}
// End namespace SpaceBallistics
