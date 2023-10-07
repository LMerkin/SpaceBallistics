// vim:ts=2:et
//===========================================================================//
//                           "MomentsOfInertia.hpp":                         //
//       Computation of MoIs for Conical and Spherical Shells, etc           //
//===========================================================================//
#pragma  once
#include "Types.hpp"
#include <cmath>
#include <cassert>
#include <cfloat>

namespace SpaceBallistics
{
  namespace Bits
  {
    //=======================================================================//
    // Pade Approximations for Sin and Cos:                                  //
    //=======================================================================//
    // We want all functions below to be "constexpr", but std::sin and std::cos
    // are not "constexpr" funcs for C++ < 26, so we provide our own, reasonably
    // accurate, implementation, using rational Pade approximations.
    // Absolute error is < 2e-10 which is OK for our applications:
    //
    //-----------------------------------------------------------------------//
    // "SinPade": Arg is assumed to be in [0..Pi/4]:                         //
    //-----------------------------------------------------------------------//
    constexpr double SinPade(double a_x)
    {
      assert(0 <= a_x && a_x < M_PI_4 + 100.0 * DBL_EPSILON);
      constexpr double a5 =  12671.0 /  4363920.0;
      constexpr double a3 = - 2363.0 /    18183.0;
      constexpr double b6 =    121.0 / 16662240.0;
      constexpr double b4 =    601.0 /   872784.0;
      constexpr double b2 =    445.0 /    12122.0;
      double           x2 = a_x * a_x;
      return a_x * ( (a5 * x2 + a3) * x2 + 1.0)
                 / (((b6 * x2 + b4) * x2 + b2) * x2 + 1.0);
    }

    //-----------------------------------------------------------------------//
    // "CosPade": Arg is assumed to be in [0..Pi/4]:                         //
    //-----------------------------------------------------------------------//
    constexpr double CosPade(double a_x)
    {
      assert(0 <= a_x && a_x < M_PI_4 + 100.0 * DBL_EPSILON);
      constexpr double a4 =  30257.0 / 1577520.0;
      constexpr double a2 = -  425.0 /     939.0;
      constexpr double b6 =     59.0 / 3155040.0;
      constexpr double b4 =   1907.0 / 1577520.0;
      constexpr double b2 =     89.0 /    1878.0;
      double           x2 = a_x * a_x;
      return ( (a4 * x2 + a2) * x2 + 1.0)  /
             (((b6 * x2 + b4) * x2 + b2) * x2 + 1.0);
    }

    //-----------------------------------------------------------------------//
    // "Sin" for an arbitrary arg:                                           //
    //-----------------------------------------------------------------------//
    constexpr double Sin(double a_x)
    {
      // First, normalise the arg to the interval [0; +00):
      bool chSgn = (a_x < 0.0);
      a_x        = std::fabs(a_x);

      // Then normalise it to the interval [0..2*Pi); NB: "fmod" is "constexpr"
      // function since C++23:
      a_x = std::fmod(a_x, 2.0 * M_PI);

      // Then to the interval [0..Pi]:
      if (a_x > M_PI)
      {
        a_x  -= M_PI;
        chSgn = !chSgn;
      }
      // Then to the interval [0..Pi/2] (same sign):
      if (a_x > M_PI_2)
        a_x = M_PI - a_x;

      // Finally, to the interval [0..Pi/4], and compute the function:
      double res = (a_x <= M_PI_4) ? SinPade(a_x) : CosPade(M_PI_2 - a_x);

      // Don't forget the sign:
      return chSgn ? (- res) : res;
    }

    //-----------------------------------------------------------------------//
    // "Cos" for an arbitrary arg:                                           //
    //-----------------------------------------------------------------------//
    constexpr double Cos(double a_x)
    {
      // First, normalise the arg to the interval [0; +00):
      a_x = std::fabs(a_x);

      // Then normalise it to the interval [0..2*Pi). NB: "fmod" is "constexpr"
      // function since C++23:
      a_x = std::fmod(a_x, 2.0 * M_PI);

      // Then to the interval [0..Pi]:
      bool chSgn = (a_x > M_PI);
      if  (chSgn)
        a_x -= M_PI;

      // Then to the interval [0..Pi/2]:
      if (a_x > M_PI_2)
      {
        a_x   = M_PI - a_x;
        chSgn = !chSgn;
      }
      // Then to the interval [0..Pi/4], and compute the function:
      double res = (a_x <= M_PI_4) ? CosPade(a_x) : SinPade(M_PI_2 - a_x);

      // Don't forget the sign:
      return chSgn ? (- res) : res;
    }

    //=======================================================================//
    // Side Surface Areas and Volumes of Construction Elements:              //
    //=======================================================================//
    // Truncated Cone:
    // "d0" and "d1" are the diameters of the bases, "h" is the height:
    //
    constexpr Area SideSurfArea_TrCone   (Len_m a_d0, Len_m a_d1, Len_m a_h)
      { return M_PI_2 * a_h * (a_d0 + a_d1); }

    constexpr Vol  Volume_TrCone         (Len_m a_d0, Len_m a_d1, Len_m a_h)
      { return (M_PI/12.0) * a_h * (Sqr(a_d0) + a_d0 * a_d1 + Sqr(a_d1)); }

    // Spherical Segment:
    // "d" is the base diameter (NOT the sphere diameter!), "h" is the height;
    // h <= d/2:
    constexpr Area SideSurfArea_SpherSegm(Len_m a_d, Len_m a_h)
    {
      assert(a_h <= a_d/2.0);
      return M_PI * (Sqr(a_d/2.0) + Sqr(a_h));
    }

    constexpr Vol Volume_SpherSegm      (Len_m a_d, Len_m a_h)
    {
      assert(a_h <= a_d/2.0);
      return M_PI * a_h * (Sqr(a_d)/8.0 + Sqr(a_h)/6.0);
    }
  }
  // End namespace "Bits"

  //=========================================================================//
  // Construction Element:                                                   //
  //=========================================================================//
  // General characteristics of standard construction elements:
  // XXX: Here we do not explicitly specify the co-ord system (typically a cer-
  // tain embedded one) in which the the Center of Mass and the Moment of Iner-
  // tia (wrt to some axis) are given. It is assumed that the co-ord system and
  // the axis are uniquely determined by the context:
  //
  struct ConstrElement
  {
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    Mass_kg const   m_mass;     // Mass
    Volume  const   m_vol;      // Notional volume (with imaginary bases)
    Len_m   const   m_CoM[3];   // Center of Mass (in some co-ord system
    MoI     const   m_MoI;      // Moment of Inertia

    //-----------------------------------------------------------------------//
    // Default Ctor:                                                         //
    //-----------------------------------------------------------------------//
    ConstrElement()
    : m_mass(0.0),
      m_vol (0.0),
      m_CoM {0.0, 0.0, 0.0},
      m_MoI (0.0)
    {}

    //-----------------------------------------------------------------------//
    // Generators:                                                           //
    //-----------------------------------------------------------------------//
    //-----------------------------------------------------------------------//
    // "TrConeXY":                                                           //
    //-----------------------------------------------------------------------//
    // Object  : Truncated (in general) conical shell (side surface only).
    // Geometry: Base diameters "d0" and "d1", height "h";  may have  d0<=>d1 ;
    //           either "d0" or "d1" may be be 0 (ie full cone), but  not both.
    // Location: The cone axis is lying in the OXY plane, at the angle  "alpha"
    //           to the positive direction of the OX axis; we can always assume
    //           that |alpha| <= Pi/2;
    //           the "d0" base diameter (NB: the resp Y-coord is irrelevant be-
    //           cause the MoI is computed  wrt the OY axis, hence not specifi-
    //           ed); the "d1" diameter corresponds to the other end of the co-
    //           ne axis segment     (where X = x1 >= x0, where x1==x0 only if
    //           |alpha|==Pi/2).
    //           When alpha=0, the truncated cone axis coincides with OX:
    // Return value: "ConstElement" corresponding to SurfDensity=1:
    //
    static ConstrElement TrConeXY
    (
      Len_m   a_x0,
      double  a_alpha, // Dimension-less, in [-Pi/2 .. Pi/2]
      Len_m   a_d0,    // Base diameter at a_x0
      Len_m   a_d1,    // Base diameter at the other section (X = x1 >= x0)
      Len_m   a_h      // Must be > 0
    )
    {
      assert(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
             !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) <= M_PI_2);
      double cosA = Cos(a_alpha);
      auto   r    = a_d0 / 2.0;
      auto   R    = a_d1 / 2.0;
      auto   R2   = Sqr(R);
      auto   r2   = Sqr(r);
      auto   th2  = 2.0 * Sqr(a_h);
      auto   s    = R + r;
      // Moment of Inertia per Surface Density. Must be strictly positive:
      auto   L4   =
        (M_PI/4.0) * SqRt(Sqr(R - r) + Sqr(a_h)) *
        (
          cosA * ((16.0/3.0) * (R + r/2.0) * a_h * a_x0 -
                  cosA * (R2 * s + (r2 - th2)  * R + (r2 - th2/3.0) * r)) +
          2.0 * s * (R2 + r2 + 2.0 * Sqr(a_x0))
        );
      assert(IsPos(L4));

      auto S = Bits::SideSurfArea_TrCone(a_d0, a_d1, a_h);
      auto V = Bits::Volume_TrCone      (a_d0, a_d1, a_h);

      return
      {
        S  * SurfDens(1.0),
        V,
        {},
        L4 * SurfDens(1.0)
      };
    }

    //-----------------------------------------------------------------------//
    // Addition:                                                             //
    //-----------------------------------------------------------------------//
    // XXX: The "right" arg is passed by copy -- for such a small struct, it is
    // probbaly better than passing it by const ref:
    MoISV operator+ (MoISV a_right) const
    {
      return MoISV
      (
        this->m_specMoI  + a_right.m_specMoI,
        this->m_S        + a_right.m_S,
        this->m_V        + a_right.m_V
      );
    }

    MoISV& operator+=(MoISV a_right)
    {
      this->m_specMoI  += a_right.m_specMoI;
      this->m_S        += a_right.m_S;
      this->m_V        += a_right.m_V;
      return *this;
    }
  };

  //=========================================================================//
  // "MoISV_TrCone_XZ":                                                      //
  //=========================================================================//
  // The geometry is similar to that in "MoISV_TrConeXY", but the axis of the
  // truncated cone is now located in the OXZ plane, at the angle "alpha" to
  // the OX axis (|alpha| <= Pi/2);  (x0,z0) are the coords of the "initial"
  // point of the axis segment (with the lowest X-coord), corresp to the base
  // diameter "d0"; "d1" is the other base's diameter. NB: in this case, both
  // "x0" and "z0" are required. The result is invariant under the transform
  // z0 <-> (-z0), alpha <-> (-alpha).
  // When z0=0 and alpha=0, we have the same case i as "MoISV_TrConeXY"  with
  // alpha=0 (the truncated cone axis coinciding with OX) (checked with Maple)!
  //
  inline MoISV MoISV_TrCone_XZ
  (
    Len_m   a_x0,          // Smallest X-coord of the base center
    Len_m   a_z0,          // Corresp  Z-coord of the base center
    double  a_alpha,       // Dimension-less, in [-Pi/2 .. Pi/2]
    Len_m   a_d0,          // Base diameter at (a_x0, a_z0)
    Len_m   a_d1,          // Base diameter at the other section (X = x1 >= x0)
    Len_m   a_h            // Must be > 0
  )
  {
    assert(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
           !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) <= M_PI_2);
    double cosA = Cos(a_alpha);
    double sinA = Sin(a_alpha);
    auto   r    = a_d0 / 2.0;
    auto   R    = a_d1 / 2.0;
    auto   R2   = Sqr(R);
    auto   r2   = Sqr(r);
    auto   s    = R + r;
    // Moment of Inertia per Surface Density. Must be strictly positive:
    auto   L4   =
      (M_PI/4.0) * SqRt(Sqr(R - r) + Sqr(a_h)) *
      (
        a_h * ((16.0/3.0)  * (R + r/2.0) * (a_x0 * cosA + a_z0 * sinA) +
               2.0 * a_h   * (R + r/3.0))                              +
        s * (R2 + r2 + 4.0 * (Sqr(a_x0) + Sqr(a_z0)))
      );
    assert(IsPos(L4));

    // Surface Area:
    auto S = SideSurfArea_TrCone(a_d0, a_d1, a_h);

    // Center of Mass:
    return MoISV(L4, S, xC, 0.0, zC);
  }

  //=========================================================================//
  // "MoISV_SpherSegm_XY":                                                    //
  //=========================================================================//
  // Similar to "MoISV_TrCone_XY", but for a Spherical Segment of the base diam-
  // eter "d" and the height (from the base plane to the pole) "h",  where  we
  // assume h <= d/2 (the equality corresponds to the case of a HemiSphere).
  // NB: this is a Shperical Segment, NOT a Spherical Slice, ie it always cont-
  // ains a pole. The Spherical Slice would be a more general case and a more
  // close analogy of the Truncated Cone, but would (arguably) have almost no
  // real applications.
  // Furthermore, "alpha" is the angle between the segment's axis and the posi-
  // tive direction of the X axis; |alpha| <= Pi/2; "x0" is the X-coord of the
  // base center (NOT of the pole; the corresp Y-coord is again irrelevant);
  // the spherical segment may be facing towards the positive direction  of the
  // OX axis (ie the X-coord of the pole is > x0), or otherise, as given by the
  // "is_pos_facing" flag.
  // When alpha=0, the axis of the Spherical Segment coincides with OX, which
  // is the most common case:
  //
  constexpr MassParams MoISV_SpherSegm_XY
  (
    bool    a_is_pos_facing,
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    Len_m   a_d,            // Base diameter
    Len_m   a_h             // Height: 0 < a_h <= a_d/2.0
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_d) && IsPos(a_h) && a_h <= a_d/2.0 &&
           std::fabs(a_alpha) <= M_PI_2);
    double cosA = Cos(a_alpha);
    double sgn  = a_is_pos_facing ? 1.0 : -1.0;
    auto   xP   = a_x0  + sgn * cosA * a_h;     // X-coord of the pole
    auto   r    = a_d / 2.0;                    // Base   radius
    auto   R    = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius
    // Therefore, R >= r >= h

    // NB: Formally, "L4" is invariant under the transform
    // isPosFacing -> !isPosFacing, alpha -> alpha + Pi, as expected;
    // however, in reality,  we always have |alpha| <= Pi/2:
    auto   L4   =
      M_PI * R * a_h *
      (2.0 * (Sqr(xP)   + a_h * (R - sgn  * xP   * cosA)) -
       a_h * ((2.0/3.0) * a_h + (R - a_h) * cosA * cosA));
    assert(IsPos(L4));

    auto S = SideSurfArea_SpherSegm(a_d, a_h);
    auto V = Volume_SpherSegm      (a_d, a_h);

    return MoISV(L4, S, V);
  }

  //=========================================================================//
  // "MoISV_SpherSegm_XZ":                                                   //
  //=========================================================================//
  // Similar to "MoISV_SpherSegm_XY", but the axis of the spherical segment is
  // in the OXZ plane, at the angle "alpha" to the positive direction  of the
  // OX axis, and the coords of the base center in that plane are (x0, z0).
  // When alpha=0 and z0=0, we have the same case as "MoISV_SpherSegm_XY" with
  // alpha=0, ie the axis of the Spherical Segment coincides with OX (checked
  // with Maple!). However, in the "_XZ" case, default args are not provided:
  //
  inline MoISV MoISV_SpherSegm_XZ
  (
    bool    a_is_pos_facing,
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    Len_m   a_z0,           // Z-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    Len_m   a_d,            // Base diameter
    Len_m   a_h             // Height: 0 < a_h <= a_d/2.0
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_h) && IsPos(a_d) && a_h <= a_d/2.0 &&
           std::fabs(a_alpha) <= M_PI_2);
    double cosA = Cos(a_alpha);
    double sinA = Sin(a_alpha);
    double sgn  = a_is_pos_facing ? 1.0 : -1.0;
    auto   xP   = a_x0  + sgn * cosA * a_h;     // X-coord of the pole
    auto   zP   = a_z0  + sgn * sinA * a_h;     // Z-coord of the pole
    auto   r    = a_d / 2.0;                    // Base   radius
    auto   R    = (Sqr(r) / a_h + a_h) / 2.0;   // Sphere radius

    // NB: Formally, "L4" is invariant under the transform
    // isPosFacing -> !isPosFacing, alpha -> alpha + Pi, as expected;
    // however, in reality,  we always have |alpha| <= Pi/2:
    auto   L4   =
      M_PI * R * a_h *
      (
        a_h * (a_h/3.0 - sgn * 2.0 * (xP * cosA + zP * sinA) + R) +
        2.0 * (Sqr(xP) + Sqr(zP))
      );
    assert(IsPos(L4));

    auto S = SideSurfArea_SpherSegm(a_d, a_h);
    auto V = Volume_SpherSegm      (a_d, a_h);

    return MoISV(L4, S, V);
  }
}
// End namespace SpaceBallistics
