// vim:ts=2:et
//===========================================================================//
//                           "MomentsOfInertia.cpp":                         //
//       Computation of MoIs for Conical and Spherical Shells, etc           //
//===========================================================================//
#include "MomentsOfInertia.h"
#include <cmath>
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "MoIS_TrCone_XY":                                                       //
  //=========================================================================//
  MoIS MoIS_TrCone_XY
  (
    Len_m   a_d0,     // Base diameter at a_x0
    Len_m   a_d1,     // Base diameter at the other section (at x = x1 >= x0)
    Len_m   a_h,      // Must be > 0
    Len_m   a_x0,
    double  a_alpha   // Dimension-less, in [-Pi/2 .. Pi/2]
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
    // Area of the truncated cone's side surface. Must also be positive:
    auto   S    = M_PI * a_h * s;
    assert(IsPos(S));
    return MoIS(L4, S);
  }

  //=========================================================================//
  // "MoIS_TrCone_XZ":                                                       //
  //=========================================================================//
  MoIS MoIS_TrCone_XZ
  (
    Len_m   a_d0,     // Base diameter at (a_x0, a_z0)
    Len_m   a_d1,     // Base diameter at the other section (at x = x1 >= x0)
    Len_m   a_h,      // Must be > 0
    Len_m   a_x0,
    Len_m   a_z0,
    double  a_alpha   // Dimension-less, in [-Pi/2 .. Pi/2]
  )
  {
    assert(IsPos(a_h)     && !IsNeg(a_d0)  && !IsNeg(a_d1) &&
           !(IsZero(a_d0) && IsZero(a_d1)) && std::fabs(a_alpha) <= M_PI_2);
    double cosA = std::cos(a_alpha);
    double sinA = std::sin(a_alpha);
    auto   r    = a_d0 / 2.0;
    auto   R    = a_d1 / 2.0;
    auto   s    = R + r;
    // Moment of Inertia per Surface Density. Must be strictly positive:
    auto   L4   =
      (M_PI/4.0) * SqRt(Sqr(R - r) + Sqr(a_h)) *
      (
        a_h * ((16.0/3.0) * (R + r/2.0) * (a_x0 * cosA + a_z0 * sinA) +
               2.0 * a_h  * (R + r/3.0))                              +
        s * (Sqr(R) + Sqr(r) + 4.0 * (Sqr(a_x0) + Sqr(a_z0)))
      );
    assert(IsPos(L4));
    // The size surface area is as in the "_XY" case:
    auto   S    = M_PI * a_h * s;
    assert(IsPos(S));
    return MoIS(L4, S);
  }

  //=========================================================================//
  // "MoIS_SpherSegm_XY":                                                    //
  //=========================================================================//
  MoIS MoIS_SpherSegm_XY
  (
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    bool    a_is_pos_facing
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
    auto   L4   =
      M_PI * R * a_h *
      (2.0 * (Sqr(xP)   + a_h * (R - sgn  * xP   * cosA)) -
       a_h * ((2.0/3.0) * a_h + (R - a_h) * cosA * cosA));
    assert(IsPos(L4));
    // The surface area of the Spherical Segment:
    auto   S    = (2.0*M_PI) * R * a_h;
    return MoIS(L4, S);
  }

  //=========================================================================//
  // "MoIS_SpherSegm_XZ":                                                    //
  //=========================================================================//
  MoIS MoIS_SpherSegm_XZ
  (
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    Len_m   a_z0,           // Z-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    bool    a_is_pos_facing
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
    auto   L4   =
      M_PI * R * a_h *
      (
        a_h * (a_h/3.0 - sgn * 2.0 * (xP * cosA + zP * sinA) + R) +
        2.0 * (Sqr(xP) + Sqr(zP))
      );
    assert(IsPos(L4));
    // The surface area of the Spherical Segment:
    auto   S    = (2.0*M_PI) * R * a_h;
    return MoIS(L4, S);
  }
}
// End namespace SpaceBallistics
