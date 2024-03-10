// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CE/ToricSegms.hpp":                  //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CE/ConstrElement.hpp"
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
      // Segm, and the Nominal Enclosed Volume:
      auto res    = PropMoICompsLow<true>(a_h);
      m_JPL0      = std::get<0>(res);
      m_JPL1      = std::get<1>(res);
      m_KPL       = std::get<2>(res);
      Vol enclVol = std::get<3>(res);
      assert(IsPos(m_JPL0) && IsPos(m_JPL1) && IsPos(m_KPL) && IsPos(enclVol));

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
    constexpr Len PropLevelOfVol(Vol a_v) const
    {
      return
        m_facingUp
        ? PropLevelOfVolUp (a_v)
        : PropLevelOfVolLow(a_v);
    }

  private:
    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    constexpr Len PropLevelOfVolLow(Vol a_v) const
    {
      assert(!IsNeg(a_v));

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

      // And finally, the Propellant level:
      return (1.0 - z) * m_R;
    }

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    // Using the invariant
    // V_up(l) + V_low(h-l) = EnclVol:
    //
    constexpr Len PropLevelOfVolUp(Vol a_v) const
    {
      assert(!IsNeg(a_v) && a_v <= GetEnclVol());
      Len res = GetHeight() - PropLevelOfVolLow(GetEnclVol() - a_v);
      assert(!IsNeg(res) && res <= GetHeight());
      return res;
    }

  public:
    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr std::tuple<Len5, Len5, Len4> PropMoIComps(Len a_l) const
    {
      assert(!IsNeg(a_l) && a_l <= GetHeight());
      return
        m_facingUp
        ? PropMoICompsUp        (a_l)
        : PropMoICompsLow<false>(a_l);
    }

    //-----------------------------------------------------------------------//
    // For the Low-Facing "ToricSegm":                                       //
    //-----------------------------------------------------------------------//
    // NB: The return type depends on the template flag. XXX: This method is
    // also called from the "ToricSegm" Ctor, when the parent classes are not
    // yet initialised -- so it MUST NOT call any methods on "RotationBody":
    //
    template<bool WithVol>
    constexpr
      std::conditional_t
      <
        WithVol,
        std::tuple<Len5, Len5, Len4, Vol>,
        std::tuple<Len5, Len5, Len4>
      >
    PropMoICompsLow(Len a_l) const
    {
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

      Len5 JP0   =
        m_pQR4 * (2.5   * acx + s * (x3 - x2/3.0 - (5.0/6.0) * x - 2.5));
      Len5 JP1   =
        m_pQR4 * (m_q2L * acx - s * cx * (m_q2L  + x - 0.5 * x2));
      Len4 KP    =
         m_pQR3 * (2.0  * acx - s * (x + 1.0) * (2.0 - (4.0/3.0) * x));

      if constexpr(WithVol)
      {
        Vol VP   = m_NV * (acx - cx  * s);
        return {JP0, JP1, KP, VP};
      }
      else
        return {JP0, JP1, KP};
    }

    //-----------------------------------------------------------------------//
    // For the Up-Facing "ToricSegm":                                        //
    //-----------------------------------------------------------------------//
    constexpr std::tuple<Len5, Len5, Len4> PropMoICompsUp(Len a_l) const
    {
      // In this case, the formulas are more involved,  so we first compute the
      // MoI components of the corresp complementary Segment, take a difference
      // wrt a Fully-Propellant-Loaded Segment, and apply a shift:
      //
      // MoI components for a complementary segment  can be obtained using the
      // same formulas as for a Low-Facing one (for which the Propellant is ad-
      // jacent to the Pole, not to the Bottom Plane):
      //
      auto [JC0, JC1, KC, VC] = PropMoICompsLow<true>(GetHeight() - a_l);

      // Subtract the above vals from the MoI components for a Fully-Loaded
      // Low-Facing Segment:
      Len5 JP0 = m_JPL0 - JC0;
      Len5 JP1 = m_JPL1 - JC1;
      Len4 KP  = m_KPL  - KC;
      Vol  VP  = GetEnclVol() - VC;
      assert(IsPos(JP0) && IsPos(JP1) && IsPos(KP) && IsPos(VP));

      // "JC0" is now almost what we need:
      // for our  required segm,  Base    @ 0,     Surface @ l,
      // for what we got so far,  Surface @ (h-l), Base    @ h;
      // the orientation does not matter, but the distance to the OEta axis is;
      // so move it in the Low direction of Xi by "h":
      JP0 += m_h2 * VP - m_twoH * KP;
      assert(IsPos(JP0));

      // "JP1" is taken as is, because the Eta co-ords are unaffected by Xi
      // shifts:
      assert(IsPos(JP1));

      // "KP" is subject to a similar shift and sign inversion (due to mirror
      // symmetry wrt OEta axis):
      KP = GetHeight() * VP - KP;
      assert(IsPos(KP));

      // All Done:
      return {JP0, JP1, KP};
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
    Len2    m_cJP0;
    Len4    m_cJP1;
    Len2    m_cKP;

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
      m_cJP0   = (Pi<double> / 3.0) * (R2 - r2);
      m_cJP1   = Pi_4<double>       * (R4 - r4);
      m_cKP    = Pi_2<double>       * (R2 - r2);

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
    constexpr Len PropLevelOfVol(Vol a_v) const
    {
      assert(!IsNeg(a_v) && a_v <= GetEnclVol());
      return double(a_v / GetEnclVol()) * GetHeight();
    }

    //=======================================================================//
    // MoI Components for the Propellant of Given Level:                     //
    //=======================================================================//
    constexpr std::tuple<Len5, Len5, Len4> PropMoIComps(Len a_l) const
    {
      assert(!IsNeg(a_l) && a_l <= GetHeight());
      Len2 l2 = Sqr(a_l);
      Len3 l3 = a_l * l2;

      Len5 JP0 = m_cJP0 * l3;
      Len5 JP1 = m_cJP1 * a_l;
      Len4 KP  = m_cKP  * l2;

      return {JP0, JP1, KP};
    }
  };
}
// End namespace SpaceBallistics
