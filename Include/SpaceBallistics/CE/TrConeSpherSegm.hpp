// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CE/TrConeSpherSegm.hpp":               //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CE/ConstrElement.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TrCone" Class: Truncated Cone Shell (w/o and with Propellant):         //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only), with
  //           SurfDensity assumed to be 1 (this may later be changed by cal-
  //           ling "ProRateMass" on the object constructed).
  // Geometry: Base diameters "d0" and "d1", height "h";  may have  d0<=>d1 ;
  //           either "d0" or "d1" may be be 0 (ie full cone), but  not both.
  // Location: XXX: Either "z0" or "y0" must be 0. The cone axis is lying in
  //           the XY or XZ plane, resp., at the angle  "alpha" to the posit-
  //           ive direction of the OX axis; we assume that |alpha| < Pi/2;
  //           the "d0" base diameter corresponds to the base center (x0,y0)
  //           (if z0=0, ie XY plane) or (x0,z0) (if y0=0, ie XZ plane); the
  //           "d1" diameter corresponds to  the other end of the cone axis
  //           segment (with X = x1 > x0);
  //           when alpha=0, the axis of the Spherical Segment coincides with
  //           OX, which is the most common case.
  // Mass:     If the "mass" param is set, its value is used  for the element
  //           mass which is in this case final. If this param is 0, the mass
  //           is set under the assumtion of SutfDensity=1, and it is NOT fi-
  //           nal yet (needs to be set later via "ProRateMass").
  // Return value: The resulting "ConstrElement" object.
  //
  class TrCone final: public RotationBody
  {
  private:
    //=======================================================================//
    // Data Flds: Truncated Cone's Geometry:                                 //
    //=======================================================================//
    // Over-all sizes:
    Len       m_r;   // Left  (Upper, Smaller-X) base radius
    Len       m_R;   // Right (Lower, Larger-X)  base radius
    Len       m_h;   // Duplicate of RotationBody::m_h which is "private"

    // Memoise pre-computed coeffs of the Cardano formula in "LevelOfVol":
    Len       m_deltaR;
    Len3      m_clVol;
    Len2      m_Rh;
    Len6      m_Rh3;

    // Default Ctor is deleted:
    TrCone() = delete;

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr static Len LevelOfVol(Vol a_v, RotationBody const* a_ctx)
    {
      // "a_ctx" must actually be a ptr to "TrCone":
      TrCone const* trc = static_cast<TrCone const*>(a_ctx);
      assert(trc != nullptr);
      return trc->LevelOfVol(a_v);
    }

    constexpr Len LevelOfVol(Vol a_v) const
    {
      assert(!IsNeg(a_v) && a_v <= GetEnclVol());
      return
        IsZero(m_deltaR)
        ? // R==r: The simplest and the most common case: A Cylinder:
          double(a_v / GetEnclVol()) * m_h

        : // General case: Solving a cubic equation by the Cardano formula;
          // since Vol'(l) > 0 globally, there is only 1 real root:
          (m_Rh - CbRt(m_Rh3 - m_clVol * a_v)) / m_deltaR;
    }

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr TrCone
    (
      Len         a_x0,       // Left  (upper) base center (Left rot.axis end)
      Len         a_y0,       //
      Len         a_z0,       //
      double      a_alpha,    // Dimension-less, in (-Pi/2 .. Pi/2)
      Len         a_d0,       // Left  (Upper, Smaller-X) base diameter
      Len         a_d1,       // Right (Lower, Larger-X)  base diameter
      Len         a_h,        // Height
      Density     a_rho,      // Propellant Density (may be 0)
      Mass        a_empty_mass = UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // Over-All Geometry:                                                  //
      //---------------------------------------------------------------------//
      assert( IsPos(a_h)   && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
            !(IsZero(a_d0) && IsZero(a_d1)) && Abs(a_alpha) < Pi_2<double> &&
            ! IsNeg(a_rho) && !IsNeg(a_empty_mass));

      // The over-all sizes:
      m_r        = a_d0 / 2.0;   // Left  (Upper, Smaller-X) base radius
      m_R        = a_d1 / 2.0;   // Right (Lower, Larger-X)  base radius
      m_h        = a_h;
      m_deltaR   = m_R - m_r;
      Len2  h2   = Sqr(m_h);

      // Memoised coeffs of the Cardano formula (used in "LevelOfVol"):
      m_clVol    = (3.0 / Pi<double>) * h2 * m_deltaR;
      m_Rh       = m_R * m_h;
      m_Rh3      = Cube (m_Rh);

      //---------------------------------------------------------------------//
      // Parent Classes Initialisation:                                      //
      //---------------------------------------------------------------------//
      // For Optimisation:
      Len    s    = SqRt(Sqr(m_deltaR) + h2);
      double a    = double(m_deltaR /  m_h);
      double a2   = Sqr(a);
      double a3   = a2 * a;
      double a4   = Sqr(a2);
      Len2   R2   = Sqr(m_R);
      Len3   R3   = R2 * m_R;
      Len4   R4   = Sqr(R2);
      Len2   r2   = Sqr(m_r);

      // Side Surface Area and Nominal Enclosed Volume:
      Area   sideSurfArea = Pi<double>     * s   * (m_R + m_r);
      Vol    enclVol      = Pi<double>/3.0 * m_h * (R2  + m_R * m_r + r2);

      // "Intrinsic" "empty" MoIs:
      Len4   JE0  =  Pi<double> * h2 * s * (m_r / 2.0  + m_R / 6.0);
      Len4   JE1  =  Pi<double>/4.0  * s * (m_R + m_r) * (R2 + r2);
      Len3   KE   = -Pi<double>/3.0  * s *  m_h * (2.0 * m_r + m_R);

      // Coeffs of "intrtinsic" MoIs with Propellant:
      double JP05 =  Pi<double> / 5.0 * a2;
      Len    JP04 = -Pi<double> / 2.0 * a  * m_R;
      Len2   JP03 =  Pi<double> / 3.0 * R2;

      double JP15 =  Pi<double> /20.0 * a4;
      Len    JP14 = -Pi<double> / 4.0 * a3 * m_R;
      Len2   JP13 =  Pi<double> / 2.0 * a2 * R2;
      Len3   JP12 = -Pi<double> / 2.0 * a  * R3;
      Len4   JP11 =  Pi<double> / 4.0 * R4;

      double KP4  = -Pi<double> / 4.0 * a2;
      Len    KP3  =  Pi<double> * 2.0 / 3.0 * a * m_R;
      Len2   KP2  = -Pi<double> / 2.0 * R2;

      // Initialise the Parent Classes' Flds: NB: For (x0,y0,z0), IsLeft=true:
      RotationBody::Init
        (sideSurfArea, enclVol, a_empty_mass,
         a_alpha,      a_x0,    a_y0,  a_z0,  true, m_h, JE0, JE1, KE,
         a_rho,        LevelOfVol,
         JP05,         JP04,    JP03,
         JP15,         JP14,    JP13,  JP12,  JP11,
         KP4,          KP3,     KP2);
    }

    //=======================================================================//
    // Non-Default Ctor, Simple Cases (TrCone/Cylinder with OX axis):        //
    //=======================================================================//
    // "TrCone" with the rotation axis coinsiding with Ox: y0=z0=alpha=0:
    //
    constexpr TrCone
    (
      Len         a_x0,       // Left  (Upper, Smaller-X) base center
      Len         a_d0,       // Left  (Upper, Smaller-X) base diameter
      Len         a_d1,       // Right (Lower, Larger-X)  base diameter
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = UnKnownMass
    )
    : TrCone(a_x0, 0.0_m, 0.0_m, 0.0, a_d0, a_d1, a_h, a_rho, a_empty_mass)
    {}

    // As above, but with d0=d1, ie a Cylinder:
    constexpr TrCone
    (
      Len         a_x0,       // Left  (Upper, Smaller-X) base center
      Len         a_d0,       // Left  (Upper, Smaller-X) base diameter
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = UnKnownMass
    )
    : TrCone(a_x0, 0.0_m, 0.0_m, 0.0, a_d0, a_d0, a_h, a_rho, a_empty_mass)
    {}
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
  // (x0, y0, z0) are the co-ords of the base center (NOT of the pole!), and
  // either "z0" or "y0" must be 0;  in the former case, the Segment's rotat-
  // ion axis is assumed to lie in the XY plane, in the latter -- in XZ.
  // The spherical segment may be facing "right" (ie towards the positive dir-
  // ection of the OX axis), or otherise, as given by the "facing_right" flag.
  // When alpha=0, the axis of the Spherical Segment coincides with OX, which
  // is the most common case:
  //
  class SpherSegm final: public RotationBody
  {
  private:
    //=======================================================================//
    // Data Flds: Spherical Segment's Geometry:                              //
    //=======================================================================//
    // Over-all sizes and Orientation:
    Len    m_r;              // Segment Base Radius
    Len    m_R;              // Sphere  Radius
    Len3   m_R3;             // Radius^3
    bool   m_facingRight;

    // Default Ctor is deleted:
    SpherSegm() = delete;

    //=======================================================================//
    // Propellant Volume -> Propellant Level:                                //
    //=======================================================================//
    constexpr static Len LevelOfVol(Vol a_v, RotationBody const* a_ctx)
    {
      // "a_ctx" must actually be a ptr to "SpherSegm":
      SpherSegm const* segm = static_cast<SpherSegm const*>(a_ctx);
      assert(segm != nullptr);
      return
        segm->m_facingRight
        ? segm->LevelOfVolRight(a_v)
        : segm->LevelOfVolLeft (a_v);
    }

    //-----------------------------------------------------------------------//
    // For the Right-Facing "SpherSegm":                                     //
    //-----------------------------------------------------------------------//
    constexpr Len LevelOfVolRight(Vol a_v) const
    {
      // Solving the equation
      //   x^2*(3-x) = v,
      // where
      //   "x" is the Propellant level relative to "R": x=l/R, 0 <= x <= 1;
      //   "v" is the Volume/(Pi/3*R^3),                       0 <= v <= 2;
      // by using Halley's Method (to avoid dealing with complex roots in the
      // Cardano formula):
      // x=1/2 is a reasonble initial approximation (XXX: We may use "h"  for
      // a more accurate initial point, but this is not required yet):
      //
      double x   = 0.5;
      double tv  = 3.0 * double(a_v / m_R3) / Pi<double>;

      constexpr double Tol = 100.0 * Eps<double>;
      assert(0.0 <= tv && tv < 2.0 + Tol);
      tv = std::min(2.0,  tv);       // Enforce the boundary

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

      // And finally, the level:
      return x * m_R;
    }

    //-----------------------------------------------------------------------//
    // For the Left-Facing "SpherSegm":                                      //
    //-----------------------------------------------------------------------//
    // Using the invariant
    // V_left(l) + V_right(h-l) = EnclVol:
    //
    constexpr Len LevelOfVolLeft(Vol a_v) const
    {
      assert(!IsNeg(a_v) && a_v <= GetEnclVol());
      Len res = GetHeight() - LevelOfVolRight(GetEnclVol() - a_v);
      assert(!IsNeg(res) && res <= GetHeight());
      return res;
    }

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    constexpr SpherSegm
    (
      bool        a_facing_right,
      Len         a_x0,          // Base center
      Len         a_y0,          //
      Len         a_z0,          //
      double      a_alpha,       // Dimension-less, in (-Pi/2 .. Pi/2)
      Len         a_d,           // Base diameter
      Len         a_h,           // Height
      Density     a_rho,         // Propellant Density (may be 0)
      Mass        a_empty_mass = UnKnownMass
    )
    {
      //---------------------------------------------------------------------//
      // Over-All Geometry:                                                  //
      //---------------------------------------------------------------------//
      assert(IsPos(a_d)   &&  IsPos(a_h) && Abs(a_alpha) < Pi_2<double> &&
            !IsNeg(a_rho) && !IsNeg(a_empty_mass));

      m_r           = a_d / 2.0;                         // Base   radius
      assert(a_h   <= m_r * (1.0 + 10.0 * Eps<double>)); // Import. constraint!
      m_R           = (Sqr(m_r) / a_h + a_h) / 2.0;      // Sphere radius
      m_R3          = Cube(m_R);
      m_facingRight = a_facing_right;

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

      // "Intrinsic" "empty" MoIs (NB: they do not depend on "FacingRight",
      // because JE0 is wrt the rotation axis "xi"; furthermore, the Left- and
      // Right-facing segments  are mutually-symmetric wrt any axis  orthogonal
      // to "xi", hence "JE1" must also be the same; for "KE", it is just a lu-
      // cky coincidence due to f(xi)*sqrt(1+f'(xi)^2)=R=const):
      //
      Len4 JE0    = TwoPi<double> / 3.0 * m_R * h3;
      Len4 JE1    = m_R * enclVol;
      Len3 KE     = -  Pi<double> * m_R * h2;

      // Coeffs of "intrtinsic" MoIs with Propellant. Unlike the "empty" ones
      // above, these coeffs do depend on the Segment orientation:
      Len Rmh  = m_R - a_h;
      assert(!IsNeg(Rmh));
      Len TRmh = m_R + Rmh;

      double JP05 = -Pi<double> / 5.0;
      Len    JP04 = (Pi<double> / 2.0)  * (a_facing_right ? m_R : -Rmh);
      Len2   JP03 =  a_facing_right
                     ? Len2(0.0)
                     : Pi<double> / 3.0 * TRmh * a_h;

      double JP15 =  Pi<double> / 20.0;
      Len    JP14 = (Pi<double> / 4.0)  * (a_facing_right ? -m_R : Rmh);

      Len2   JP13 =  Pi<double> *
                     (a_facing_right
                      ? R2 / 3.0
                      : R2 / 3.0 - m_R * a_h + h2 / 2.0);

      Len3   JP12 =  a_facing_right
                     ? Len3(0.0)
                     : (-Pi<double> / 2.0) * Rmh * TRmh * a_h;

      Len4   JP11 =  a_facing_right
                     ? Len4(0.0)
                     : ( Pi<double> / 4.0) * Sqr(TRmh)  * h2;

      double KP4  =  Pi<double> / 4.0;
      Len    KP3  = (Pi<double> * 2.0/3.0) * (a_facing_right  ? -m_R : Rmh);
      Len2   KP2  =  a_facing_right
                     ? Len2(0.0) : (-Pi<double> / 2.0) * TRmh * a_h;

      // Initialise the Parent Classes' Flds:
      // NB: (x0,y0,z0) is the Base Center, so for it, IsLeft = FacingRight:
      RotationBody::Init
        (sideSurfArea, enclVol, a_empty_mass,
         a_alpha,   a_x0, a_y0, a_z0, a_facing_right, a_h, JE0, JE1, KE,
         a_rho,     LevelOfVol,
         JP05,      JP04, JP03,
         JP15,      JP14, JP13, JP12, JP11,
         KP4,       KP3,  KP2);
    }

    //=======================================================================//
    // Non-Default Ctor, Simple Cases:
    //=======================================================================//
    // Rotation axis coinciding with OX:
    //
    constexpr SpherSegm
    (
      bool        a_facing_right,
      Len         a_x0,       // Base center
      Len         a_d,        // Base diameter
      Len         a_h,        // Height
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = UnKnownMass
    )
    : SpherSegm(a_facing_right, a_x0, 0.0_m, 0.0_m, 0.0, a_d, a_h, a_rho,
                a_empty_mass)
    {}

    // As above, but with d/2 = h, ie a HemiSpehere:
    //
    constexpr SpherSegm
    (
      bool        a_facing_right,
      Len         a_x0,       // Base center
      Len         a_d,        // Base diameter
      Density     a_rho,      // 0 may be OK
      Mass        a_empty_mass = UnKnownMass
    )
    : SpherSegm(a_facing_right, a_x0, 0.0_m, 0.0_m, 0.0, a_d, a_d / 2.0, a_rho,
                a_empty_mass)
    {}
  };
}
// End namespace SpaceBallistics
