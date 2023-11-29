// vim:ts=2:et
//===========================================================================//
//                             "ConstrElement.hpp":                          //
//                 Geometrical Objects as Construction Elements              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <cmath>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // Construction Element:                                                   //
  //=========================================================================//
  // General characteristics of standard construction elements:
  // Currently, 2D Shells (Truncated Cones and Spherical Segments), 3D Volumes
  // (of the same shape as Shells) and Point Masses are supported.  For them, 
  // we can compute Centers of Masses (CoM) and Moments of Inertia (MoI) wrt
  // the OX, OY and OZ axes.
  // XXX:
  // Here we do not explicitly specify the co-ord system OXYZ (typically a cer-
  // tain embedded one) in which the CoM  and the MoI (wrt to the OY axis)  are
  // given. It is assumed that the co-ord system OXYZ is uniquely determined by
  // the context:
  //
  class ConstrElement
  {
  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    Mass  m_mass;        // Mass
    Area  m_surfArea;    // Side Surface Area (w/o Bases)
    Vol   m_vol;         // Notional volume (with imaginary Bases added)
    Len   m_CoM [3];     // (X,Y,Z) co-ords of the Center of Masses
    MoI   m_MoIs[3];     // Moments of Inertia wrt the OX, OY and OZ axes
    bool  m_massIsFinal; // If not set, "m_mass" and "m_MoIs" are not valid yet

  public:
    //=======================================================================//
    // Default Ctor:                                                         //
    //=======================================================================//
    // NB: All data fields are initialised to 0s, not to NaNs, so the "empty"
    // obj (constructed by the Default Ctor) can be used as an initial value
    // for summation ("+", "+="):
    //
    constexpr ConstrElement()
    : m_mass       (0.0),
      m_surfArea   (0.0),
      m_vol        (0.0),
      m_CoM        {0.0_m,    0.0_m,    0.0_m},
      m_MoIs       {MoI(0.0), MoI(0.0), MoI(0.0)},
      m_massIsFinal(false)
    {}

    //=======================================================================//
    // Copy Ctor, Assignment, Equality:                                      //
    //=======================================================================//
    // They are trivial. Also, the Dtor is trivial, and is auto-generated:
    //
    constexpr ConstrElement            (ConstrElement const& a_right) = default;
    constexpr ConstrElement& operator= (ConstrElement const& a_right) = default;

    constexpr bool           operator==(ConstrElement const& a_right) const
      = default;
    constexpr bool           operator!=(ConstrElement const& a_right) const
      = default;

  private:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    // For internal use only. NB: The "a_mass" param may be final or provisio-
    // nal, and the "a_mass_is_valid" is set accordingly.
    // External callers should use various Generators (see below) instead:
    //
    constexpr ConstrElement
    (
      Mass a_mass,  Area a_surf_area, Vol a_vol,
      Len  a_xc,    Len  a_yc,        Len a_zc,
      MoI  a_moi_x, MoI  a_moi_y,     MoI a_moi_z,
      bool a_mass_is_final
    )
    : m_mass        (a_mass),
      m_surfArea    (a_surf_area),
      m_vol         (a_vol),
      m_CoM         {a_xc,    a_yc,    a_zc},
      m_MoIs        {a_moi_x, a_moi_y, a_moi_z},
      m_massIsFinal (a_mass_is_final)
    {
      assert(!(IsNeg(m_mass)    || IsNeg(m_surfArea) || IsNeg(m_vol)    ||
               IsNeg(m_MoIs[0]) || IsNeg(m_MoIs[1])  || IsNeg(m_MoIs[2])));
    }

  public:
    //=======================================================================//
    // "SetMass":                                                            //
    //=======================================================================//
    // XXX: A "setter" for the Mass (invoked when we know the final mass, sets
    // it and pro-rates the mass-dependent params). After that, the mass  will
    // be final and cannot be changed anymore:
    //
    constexpr void SetMass(Mass a_final_mass)
    {
      // Cannot adjust the Mass once it has been finalised:
      assert(!m_massIsFinal);

      // Bot the curr and the new mass must be STRICTLY positive (if the former
      // is not, the scaling factor would be infinite; if the latter is not, it
      // just makes no sense):
      assert(IsPos(m_mass) && IsPos(a_final_mass));

      // If all is OK:
      double scaleCoeff = double(a_final_mass / m_mass);
      assert(scaleCoeff > 0.0);

      // Set the Mass and adjust the MoIs:
      m_mass        = a_final_mass;
      m_MoIs[0]    *= scaleCoeff;
      m_MoIs[1]    *= scaleCoeff;
      m_MoIs[2]    *= scaleCoeff;
      m_massIsFinal = true;
    }

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    // Geometrical params are always available, even if the Mass is not final
    // yet:
    constexpr Area GetSurfArea() const { return m_surfArea; }
    constexpr Vol  GetVol()      const { return m_vol;      }

    using     Point = decltype(m_CoM);
    using     MoIs  = decltype(m_MoIs);

    constexpr Point const& GetCoM() const { return m_CoM;   }

    // But the Mass and MoIs  may or may not be available.
    // If not, assert failure will be signaled (NB: in C++ >= 14, "assert" is
    // allowed in "constexpr" functions, whereas throwing exceptions is not):
    //
    constexpr Mass GetMass() const
    {
      assert(m_massIsFinal);
      return m_mass;
    }

    constexpr MoIs const& GetMoIs() const
    {
      assert(m_massIsFinal);
      return m_MoIs;
    }

    //=======================================================================//
    // "ShellTrCone" Generator:                                              //
    //=======================================================================//
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
    //           nal yet (needs to be set later via "SetMass").
    // Return value: The resulting "ConstrElement" object.
    //
    constexpr static ConstrElement ShellTrCone
    (
      Len     a_x0,
      Len     a_y0,
      Len     a_z0,
      double  a_alpha,     // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d0,        // Base diameter at x0
      Len     a_d1,        // Base diameter at the other section (x1 > x0)
      Len     a_h,         // Must be > 0
      Mass    a_mass       // 0 if the mass is not known yet
    )
    {
      assert( IsPos(a_h)   && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
            !(IsZero(a_d0) && IsZero(a_d1)) && Abs(a_alpha) < M_PI_2);

      // At least one of "z0", "y0" must be 0, because the TrCone axis lies
      // either in the OXY plane or in the OXZ plane (or indeed in both, if
      // that axis coincides with OX):
      //
      bool inXY = IsZero(a_z0);
      bool inXZ = IsZero(a_y0);
      assert(inXY || inXZ);
      Len  yz0  = inXY ? a_y0 : a_z0;

      double cosA  = Cos(a_alpha);
      double sinA  = Sin(a_alpha);
      double tanA  = sinA  / cosA;

      Len    r     = a_d0  / 2.0;   // "Left"  (smaller "x") base radius
      Len    R     = a_d1  / 2.0;   // "Right" (larger  "x") base radius
      Len    s     = R + r;
      Len    s1    = R + r / 2.0;
      auto   R2    = Sqr(R);
      auto   r2    = Sqr(r);
      auto   h2    = Sqr(a_h);
      auto   th2   = 2.0 * h2;
      auto   yz02  = Sqr(yz0);

      // The coeff which is present in all MoI exprs below:
      auto   coeff = (M_PI/4.0) * SqRt(Sqr(R - r) + Sqr(a_h));

      // The common term present in "L4x" and "L4in":
      auto   cTerm = R2 * s + (r2 - th2) * R + (r2 - th2/3.0) * r;

      // Relative MoIs (per Surface Density): Denoted "L4{x,in,orth}", since
      // their dim is Len^4. Must all be strictly > 0:
      //
      // "L4x": L4 wrt OX:
      auto L4x =
        coeff *
        (cTerm * (1.0 + Sqr(cosA)) + (16.0/3.0) * yz0 * a_h * s1 * sinA +
         4.0   * (R * (h2 + yz02)  + r * (yz02  + h2/3.0)));

      // "L4in":
      // L4 wrt the other principal axis lying   IN the plane  (OX,TrConeAxis):
      // that is, if the TrConeAxis lies in OXY, it is OY;  if in OXZ, then XZ;
      // if in both (ie TrConeAxis coincides with OX), then "L4in" and "L4orth"
      // (see below) must be equal.
      // NB: "L4in" is invariant of the sign of "alpha":
      auto L4in =
        coeff *
        (cosA * ((16.0/3.0) * a_h *  s1 * a_x0 -
                cosA * (R2  * s   + (r2 - th2) * R + (r2 - th2/3.0) * r)) +
         2.0 * s * (R2 + r2 + 2.0 * Sqr(a_x0)));
      assert(IsPos(L4in));

      // "L4orth":
      // L4 wrt the axis orthogonal to the (OX,TrConeAxis) plane;  that is, if
      // the TrConeAxis lies in OXY, it is OZ; if in OXZ, then XY; if  in both
      // (ie TrConeAxis coincides with OX),  then "L4in"a and "L4orth" must be
      // equal.
      // NB: "L4orth" is invariant iunder the transform
      //              yz0 <-> (-yz0), alpha <-> (-alpha),
      // where yz0 is "a_y0" or "a_z0" depending on the plane of TrConeAxis:
      //
      auto L4orth =
        coeff *
        (a_h  * ((16.0/3.0) * s1 * (a_x0 * cosA + yz0 * sinA) +
                 2.0 * a_h  * (R + r/3.0))                    +
         s * (R2 + r2 + 4.0 * (Sqr(a_x0)  + Sqr(yz0))));
      assert(IsPos(L4orth));

      // And indeed, if alpha==0 and yz0==0, then indeed L4in == L4orth (verif-
      // ied with Maple).

      // Side Surface Area and Notional Volume:
      auto S = SideSurfArea_TrCone(a_d0, a_d1, a_h);
      auto V = Volume_TrCone      (a_d0, a_d1, a_h);

      // The mass may or may not be given. If it is given,  calculate the Surf-
      // ace Density; otherwise, assume the SurfDens to be 1.0.   Then set the
      // Mass and MoIs:
      bool hasMass  = IsPos(a_mass);
      auto surfDens = hasMass ? (a_mass / S) : SurfDens(1.0);
      Mass M        = hasMass ?  a_mass      : (S * surfDens);

      // DON'T DO IT:
      //decltype(1_kg / Sqr(1_m)) surfDens;
      //if (hasMass)
      //  surfDens = a_mass / S;
      //else
      //  surfDens = SurfDens(1.0);

      // So the actial MoIs are:
      MoI  MoIX     = L4x                      * surfDens;
      MoI  MoIY     = (inXY ? L4in   : L4orth) * surfDens;
      MoI  MoIZ     = (inXY ? L4orth : L4in)   * surfDens;

      // Co-ords of the CoM:
      // The X-coord is such that the above L4 is minimal. Other co-ords are
      // computed from the rotation axis geometry:
      Len  deltaXC  = cosA * (2.0 * R + r) / (3.0 * (R + r)) * a_h;
      Len  xC       = a_x0 + deltaXC;
      Len  yC       = inXY ? (a_y0 + deltaXC * tanA) : 0.0_m;
      Len  zC       = inXZ ? (a_z0 + deltaXC * tanA) : 0.0_m;

      return ConstrElement(M, S, V, xC, yC, zC, MoIX, MoIY, MoIZ, hasMass);
    }

    //=======================================================================//
    // "ShellSpherSegm":                                                     //
    //=======================================================================//
    // Similar to "ShellTrCone", but for a Spherical Segment of the base diam-
    // eter "d" and the height (from the base plane to the pole) "h", where we
    // assume h <= d/2 (the equality corresponds to the case of a HemiSphere).
    // NB: this is a Shperical Segment,  NOT a Spherical Slice,  ie it always
    // contains a pole. The Spherical Slice would be a more general case and a
    // more close analogy of the Truncated Cone, but it would (arguably) have
    // almost no real applications.
    // Furthermore, "alpha" is the angle between the segment's axis and the po-
    // sitive direction of the X axis; we assume that |alpha| < Pi/2;
    // (x0, y0, z0) are the co-ords of the base center (NOT of the pole!), and
    // either "z0" or "y0" must be 0;  in the former case, the Segment's rotat-
    // ion axis is assumed to lie in the XY plane, in the latter -- in XZ.
    // The spherical segment may be facing towards the positive direction   of
    // the OX axis (ie the X-coord of the pole is > x0), or otherise, as given
    // by the "is_pos_facing" flag.
    // When alpha=0, the axis of the Spherical Segment coincides with OX, which
    // is the most common case:
    //
    constexpr static ConstrElement ShellSpherSegm
    (
      bool    a_is_pos_facing,
      Len     a_x0,           // X-coord of the base center (NOT of the pole!)
      Len     a_y0,
      Len     a_z0,
      double  a_alpha,        // Dimension-less, in (-Pi/2 .. Pi/2)
      Len     a_d,            // Base diameter
      Len     a_h,            // Height: 0 < a_h <= a_d/2.0
      Mass    a_mass          // 0 if the mass is not known yet
    )
    {
      // NB: We must always have 0 < h <= R:
      assert(IsPos(a_d) && IsPos(a_h) && a_h <= a_d/2.0 &&
             Abs(a_alpha) < M_PI_2);

      // The mass may or may not be given:
      bool hasMass = IsPos(a_mass);

      // At least one of "z0", "y0" must be 0:
      bool inXY    = IsZero(a_z0);
      bool inXZ    = IsZero(a_y0);
      assert(inXY || inXZ);

      double cosA  = Cos(a_alpha);
      double sinA  = Sin(a_alpha);
      double tanA  = sinA / cosA;

      // X- and Z-coords of the Pole; the Y-coord is not used:
      double sgn   = a_is_pos_facing ? 1.0 : -1.0;
      auto   xP    = a_x0 +  sgn  * cosA * a_h;
      auto   zP    = inXZ ? (a_z0 + sgn  * sinA * a_h) : 0.0_m;
      auto   r     = a_d / 2.0;                    // Base   radius
      auto   R     = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius
      // Therefore, R >= r >= h, but we cannot formally assert this due to pos-
      // sible rounding effects...

      // NB: Formally, "L4" is invariant under the transform
      // isPosFacing -> !isPosFacing,   alpha -> alpha + Pi,
      // as one can expect; however, in reality, we always have |alpha| < Pi/2:
      auto   L4    =
        M_PI * R * a_h *
        (inXY
         ? (2.0 * (Sqr(xP)   + a_h * (R - sgn  * xP   * cosA)) -
            a_h * ((2.0/3.0) * a_h + (R - a_h) * cosA * cosA))
         : (2.0 * (Sqr(xP)   + Sqr(zP))                        +
            a_h * (a_h/3.0   - sgn * 2.0 * (xP * cosA + zP * sinA) + R))
        );
      assert(IsPos(L4));

      // Side Surface Area and Notional Volume:
      auto S       = SideSurfArea_SpherSegm(a_d, a_h);
      auto V       = Volume_SpherSegm      (a_d, a_h);

      // If the mass is given, calculate the Surface Density; otherwise, assume
      // it is 1.0. Then set the Mass and MoI:
      auto surfDens = hasMass ? (a_mass / S) : SurfDens(1.0);
      auto M        = hasMass ?  a_mass      : (S * surfDens);
      auto MoIY     = L4 * surfDens;
      auto MoIX     = MoIY;
      auto MoIZ     = MoIY;

      // Co-ords of the CoM:
      // The X-coord is such that the above L4 is minimal. Other co-ords are
      // computed from the rotation axis geometry:
      // PosFacing (sgn = 1): xC < xP;
      // NegFacing (sgn =-1): xC > xP:
      auto deltaXC = -sgn * a_h * cosA / 2.0;
      auto xC      = a_x0 + deltaXC;
      auto yC      = inXY ? (xP + deltaXC * tanA) : 0.0_m;
      auto zC      = inXZ ? (zP + deltaXC * tanA) : 0.0_m;
      assert(( a_is_pos_facing && xC < xP) ||
             (!a_is_pos_facing && xC > xP));

      return ConstrElement(M, S, V, xC, yC, zC, MoIX, MoIY, MoIZ, hasMass);
    }

    //=======================================================================//
    // "PointMass":                                                          //
    //=======================================================================//
    // A positive mass concentrated in the (x0, y0, z0) point:
    //
    constexpr static ConstrElement PointMass
      (Len a_x0, Len a_y0, Len a_z0, Mass a_mass)
    {
      assert(IsPos(a_mass));

      auto x2   = Sqr(a_x0);
      auto y2   = Sqr(a_y0);
      auto z2   = Sqr(a_z0);
      auto MoIX = a_mass * (y2 + z2); // Distance^2 to the OX axis
      auto MoIY = a_mass * (x2 + z2); // Distance^2 to the OY axis
      auto MoIZ = a_mass * (x2 + y2); // Distance^2 to the OZ axis

      // The mass is considered to be final:
      return ConstrElement
             (a_mass, Area(0.0), Vol(0.0), a_x0, a_y0, a_z0, MoIX, MoIY, MoIZ,
              true);
    }

    //=======================================================================//
    // Addition:                                                             //
    //=======================================================================//
    // The summands must both have final or non-final masses. Only one summand
    // is allowed to have zero mass, otherwise we cannot compute the CoM:
    //
    constexpr ConstrElement& operator+= (ConstrElement const& a_right)
    {
      assert(m_massIsFinal == a_right.m_massIsFinal    &&
             !(IsZero(m_mass) && IsZero(a_right.m_mass)));

      // Masses, SurfAreas, Vols and MoIs are directly-additive:
      auto m0     = m_mass;
      m_mass     += a_right.m_mass;
      m_surfArea += a_right.m_surfArea;
      m_vol      += a_right.m_vol;
      m_MoIs[0]  += a_right.m_MoIs[0];
      m_MoIs[1]  += a_right.m_MoIs[1];
      m_MoIs[2]  += a_right.m_MoIs[2];

      // For the CoM, do the weighted avg (but the total mass must be non-0):
      assert(IsPos(m_mass));
      double mu0 = double(m0             / m_mass);
      double mu1 = double(a_right.m_mass / m_mass);
      m_CoM[0]   = mu0 * m_CoM[0]  + mu1 * a_right.m_CoM[0];
      m_CoM[1]   = mu0 * m_CoM[1]  + mu1 * a_right.m_CoM[1];
      m_CoM[2]   = mu0 * m_CoM[2]  + mu1 * a_right.m_CoM[2];
      return *this;
    }

    constexpr ConstrElement operator+ (ConstrElement const& a_right) const
    {
      ConstrElement res = *this;
      res += a_right;
      return res;
    }

    //=======================================================================//
    // UTILS: Side Surface Areas and Volumes of Construction Elements:       //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Shell -- Truncated Cone:                                              //
    //-----------------------------------------------------------------------//
    // "d0" and "d1" are the diameters of the bases, "h" is the height:
    //
    constexpr static Area SideSurfArea_TrCone(Len a_d0, Len a_d1, Len a_h)
    {
      assert(!(IsNeg(a_d0) || IsNeg(a_d1) || (IsZero(a_d0) && IsZero(a_d1))) &&
             IsPos(a_h));
      return M_PI_2 * a_h * (a_d0 + a_d1);
    }

    constexpr static Vol  Volume_TrCone      (Len a_d0, Len a_d1, Len a_h)
    {
      assert(!(IsNeg(a_d0) || IsNeg(a_d1) || (IsZero(a_d0) && IsZero(a_d1))) &&
             IsPos(a_h));
      return (M_PI/12.0) * a_h * (Sqr(a_d0) + a_d0 * a_d1 + Sqr(a_d1));
    }

    //-----------------------------------------------------------------------//
    // Shell -- Spherical Segment:                                           //
    //-----------------------------------------------------------------------//
    // "d" is the base diameter (NOT the sphere diameter!), "h" is the height;
    // h <= d/2:
    constexpr static Area SideSurfArea_SpherSegm(Len a_d, Len a_h)
    {
      assert(IsPos(a_h) && a_h <= a_d/2.0);
      return M_PI * (Sqr(a_d/2.0) + Sqr(a_h));
    }

    constexpr static Vol  Volume_SpherSegm      (Len a_d, Len a_h)
    {
      assert(IsPos(a_h) && a_h <= a_d/2.0);
      return M_PI * a_h * (Sqr(a_d)/8.0 + Sqr(a_h)/6.0);
    }
  };
}
// End namespace SpaceBallistics
