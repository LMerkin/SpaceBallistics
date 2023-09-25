// vim:ts=2:et
//===========================================================================//
//                           "MomentsOfInertia.hpp":                         //
//       Computation of MoIs for Conical and Spherical Shells, etc           //
//===========================================================================//
#pragma  once
#include "Types.hpp"
#include <cmath>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // Extra Dim Types:                                                        //
  //=========================================================================//
  using Len4_m4 = decltype(IPow<4>(Len_m_1)); // MoI per Surface Density (L^4)

  //=========================================================================//
  // "MoISV" Struct:                                                         //
  //=========================================================================//
  // For convenience of using the functions below, define addition of the foll-
  // owing pairs:
  struct MoISV
  {
    Len4_m4   m_specMoI;    // MoI per Surface Density
    Area      m_S;          // Side Surface Area
    Vol       m_V;          // Nominal Volume of the Shell

    // Default Ctor sets all components to 0:
    MoISV()
    : m_specMoI(0.0),
      m_S      (0.0),
      m_V      (0.0)
    {}

    // Non-Default Ctor:
    // All components are assumed to be non-0:
    MoISV(Len4_m4 a_spec_moi, Area a_s, Vol a_v)
    : m_specMoI (a_spec_moi),
      m_S       (a_s),
      m_V       (a_v)
    { assert(IsPos(m_specMoI) && IsPos(m_S) && IsPos(m_V)); }

    // Addition:
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
  // "MoISV_TrCone_XY":                                                      //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only).
  // Geometry: Base diameters "d0" and "d1", height "h";  may have  d0<=>d1 ;
  //           either "d0" or "d1" may be be 0 (ie full cone), but  not both.
  // Location: The cone axis is lying in the OXY plane, at the angle  "alpha"
  //           to the positive direction of the OX axis; we can always assume
  //           that |alpha| <= Pi/2;
  //           "x0" is the lowest X-coord of the cone axis segment, corresp to
  //           the "d0" base diameter (NB: the resp Y-coord is irrelevant since
  //           the MoI is computed  wrt the OY axis, hence not specified); the
  //           "d1" diameter corresponds to the other end of the cone axis seg-
  //           ment (with X = x1 >= x0,  where x1==x0 only if |alpha|==Pi/2).
  //           When alpha=0, the truncated cone axis coincides with OX.
  // Result  : {Moment of Inertia wrt the OY axis per Surface Density [L^4],
  //            Side Surface Area [L^2]}:
  //
  inline MoISV MoISV_TrCone_XY
  (
    Len_m   a_d0,          // Base diameter at a_x0
    Len_m   a_d1,          // Base diameter at the other section (X = x1 >= x0)
    Len_m   a_h,           // Must be > 0
    Len_m   a_x0,
    double  a_alpha = 0.0  // Dimension-less, in [-Pi/2 .. Pi/2]
  )
  {
    assert(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
           !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) <= M_PI_2);
    double cosA = std::cos(a_alpha);
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

    // Area of the truncated cone's side surface:
    auto   S    = M_PI * a_h * s;
    // Volume of the truncated cone:
    auto   V    = (M_PI/3.0) * a_h * (R2 + R*r + r2);

    return MoISV(L4, S, V);
  }

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
  // However, in the "_XZ" case, the default args are not provided:
  //
  inline MoISV MoISV_TrCone_XZ
  (
    Len_m   a_d0,          // Base diameter at (a_x0, a_z0)
    Len_m   a_d1,          // Base diameter at the other section (X = x1 >= x0)
    Len_m   a_h,           // Must be > 0
    Len_m   a_x0,
    Len_m   a_z0,
    double  a_alpha        // Dimension-less, in [-Pi/2 .. Pi/2]
  )
  {
    assert(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
           !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) <= M_PI_2);
    double cosA = std::cos(a_alpha);
    double sinA = std::sin(a_alpha);
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

    // The side surface area and the volume are as in the "_XY" case:
    auto   S    = M_PI * a_h * s;
    auto   V    = (M_PI/3.0) * a_h * (R2 + R*r + r2);

    return MoISV(L4, S, V);
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
  inline MoISV MoISV_SpherSegm_XY
  (
    bool    a_is_pos_facing,
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    double  a_alpha = 0.0   // Dimension-less, in [-Pi/2 .. Pi/2]
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_d) && IsPos(a_h) && a_h <= a_d/2.0 &&
           std::fabs(a_alpha) <= M_PI_2);
    double cosA = std::cos(a_alpha);
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

    // The surface area and the volume of the Spherical Segment:
    auto   S    = (2.0*M_PI) * R  * a_h;
    auto   V    = M_PI * Sqr(a_h) * (R - a_h/3.0);
    assert(IsPos(V));
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
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    Len_m   a_z0,           // Z-coord of the base center (NOT of the pole!)
    double  a_alpha         // Dimension-less, in [-Pi/2 .. Pi/2]
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_h) && IsPos(a_d) && a_h <= a_d/2.0 &&
           std::fabs(a_alpha) <= M_PI_2);
    double cosA = std::cos(a_alpha);
    double sinA = std::sin(a_alpha);
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

    // The surface area and the volume of the Spherical Segment are as in the
    // "_XY" case:
    auto   S    = (2.0*M_PI) * R  * a_h;
    auto   V    = M_PI * Sqr(a_h) * (R - a_h/3.0);
    assert(IsPos(V));
    return MoISV(L4, S, V);
  }
}
// End namespace SpaceBallistics
