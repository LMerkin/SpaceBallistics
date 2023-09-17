// vim:ts=2:et
//===========================================================================//
//                           "MomentsOfInertia.hpp":                         //
//       Computation of MoIs for Conical and Spherical Shells, etc           //
//===========================================================================//
#pragma  once
#include "Types.hpp"
#include <cmath>

namespace SpaceBallistics
{
  //=========================================================================//
  // Types:                                                                  //
  //=========================================================================//
  using Len4_m4 = decltype(IPow<4>(Len_m_1)); // MoI per Surface Density (L^4)
  using Area_m2 = decltype(Sqr    (Len_m_1)); // Area (L^2)

  //-------------------------------------------------------------------------//
  // "MoIS" Struct:                                                          //
  //-------------------------------------------------------------------------//
  // For convenience of using the functions below, define addition of the foll-
  // owing pairs:
  struct MoIS: public std::pair<Len4_m4, Area_m2>
  {
    // Default Ctor sets both components to 0:
    MoIS()
    : std::pair<Len4_m4, Area_m2>(0.0, 0.0)
    {}

    // Non-Default Ctor just inherits from that of "std::pair":
    MoIS(Len4_m4 a_first, Area_m2  a_second)
    : std::pair<Len4_m4,  Area_m2>(a_first, a_second)
    {}

    // Addition:
    MoIS  operator+ (MoIS a_right) const
    {
      return MoIS
      (
        this->first  + a_right.first,
        this->second + a_right.second
      );
    }

    MoIS& operator+=(MoIS a_right)
    {
      this->first  += a_right.first;
      this->second += a_right.second;
      return *this;
    }
  };

  //=========================================================================//
  // "MoIS_TrCone_XY":                                                       //
  //=========================================================================//
  // Object  : Truncated (in general) conical shell (side surface only).
  // Geometry: Radii "R" and "r", height "h"; may have r<=>R, and either "r"
  //           or "R" may be be 0 (full cone), but not both.
  // Location: Axis in the OXY plane, at the angle "alpha"  to the OX axis;
  //           "x0" is the lowest X-coord of the axis range (the corresp Y-
  //           coord is irrelevant since the MoI is computed wrt the OY axis),
  //           corresp to the "r" radius; the "R" radius corresponds to the
  //           other end of the cone axis (with x > x0).
  // Result  : {Moment of Inertia wrt the OY axis per Surface Density [L^4],
  //            Side Surface Area [L^2]}:
  //
  inline MoIS MoIS_TrCone_XY
  (
    Len_m   a_r,
    Len_m   a_R,
    Len_m   a_h,
    Len_m   a_x0,
    double  a_alpha    // Dimension-less
  )
  {
    assert(IsPos(a_h)    && !IsNeg(a_r) && !IsNeg(a_R) &&
           !(IsZero(a_r) && IsZero(a_R)));
    double cosA = std::cos(a_alpha);
    auto   R2   = Sqr(a_R);
    auto   r2   = Sqr(a_r);
    auto   th2  = 2.0 * Sqr(a_h);
    auto   s    = a_R + a_r;
    // Moment of Inertia per Surface Density. Must be strictly positive:
    auto   L4   =
      (M_PI/4.0) * SqRt(Sqr(a_R - a_r) + Sqr(a_h)) *
      (
        cosA * ((16.0/3.0) * (a_R + a_r/2.0) * a_h * a_x0 -
                cosA * (R2 * s + (r2 - th2)  * a_R + (r2 - th2/3.0) * a_r)) +
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
  // Similar to "MoIS_TrConeXY", but the axis of the truncated cone is now loc-
  // ated in the OXZ plane, at the angle "alpha" to the OX axis, with the init-
  // ial axis point coords "x0" and "z0" (both are required in this case):
  //
  inline MoIS MoIS_TrCone_XZ
  (
    Len_m   a_r,
    Len_m   a_R,
    Len_m   a_h,
    Len_m   a_x0,
    Len_m   a_z0,
    double  a_alpha     // Dimension-less
  )
  {
    assert(IsPos(a_h)    && !IsNeg(a_r) && !IsNeg(a_R) &&
           !(IsZero(a_r) && IsZero(a_R)));
    double cosA = std::cos(a_alpha);
    double sinA = std::sin(a_alpha);
    auto   s    = a_R + a_r;
    // Moment of Inertia per Surface Density. Must be strictly positive:
    auto   L4   =
      (M_PI/4.0) * SqRt(Sqr(a_R - a_r) + Sqr(a_h)) *
      (
        a_h * ((16.0/3.0) * (a_R + a_r/2.0) * (a_x0 * cosA + a_z0 * sinA) +
               2.0 * a_h  * (a_R + a_r/3.0))                              +
        s * (Sqr(a_R) + Sqr(a_r) + 4.0 * (Sqr(a_x0) + Sqr(a_z0)))
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
  // Similar to "MoIS_TrCone_XY", but for a Spherical Segment of radius "R" and
  // height (from the section plane to the pole) "h";   "x0"  is the X-coord of
  // the pole, "alpha" is the angle of the segment's rotation axis to OX in the
  // OXY plane:
  //
  inline MoIS MoIS_SpherSegm_XY
  (
    Len_m   a_R,
    Len_m   a_h,
    Len_m   a_x0,
    double  a_alpha    // Dimension-less
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_h) && IsPos(a_R) && a_h <= a_R);
    double cosA = std::cos(a_alpha);
    auto   L4   =
      M_PI * a_R * a_h *
      (2.0 * (Sqr(a_x0) + a_h * (a_R + a_x0 * cosA)) -
       a_h * ((2.0/3.0) * a_h + (a_R - a_h) * cosA * cosA));
    assert(IsPos(L4));
    // The surface area of the Spherical Segment:
    auto   S    = (2.0*M_PI) * a_R * a_h;
    return MoIS(L4, S);
  }

  //=========================================================================//
  // "MoIS_SpherSegm_XZ":                                                    //
  //=========================================================================//
  // Similar to "MoIS_SpherSegm_XY", but the axis of the spherical segment is
  // in the OXZ plane, at the angle "alpha" to the OX axis,  and the pole co-
  // ords in that plane are (x0, z0):
  //
  inline MoIS MoIS_SpherSegm_XZ
  (
    Len_m   a_R,
    Len_m   a_h,
    Len_m   a_x0,
    Len_m   a_z0,
    double  a_alpha    // Dimension-less
  )
  {
    // NB: We must always have 0 < h <= R:
    assert(IsPos(a_h) && IsPos(a_R) && a_h <= a_R);
    double cosA = std::cos(a_alpha);
    double sinA = std::sin(a_alpha);
    auto   L4   =
      M_PI * a_R * a_h *
      (
        a_h * (a_h/3.0   + 2.0 * (a_x0 * cosA + a_z0 * sinA) + a_R) +
        2.0 * (Sqr(a_x0) + Sqr(a_z0))
      );
    assert(IsPos(L4));
    // The surface area of the Spherical Segment:
    auto   S    = (2.0*M_PI) * a_R * a_h;
    return MoIS(L4, S);
  }
}
// End namespace SpaceBallistics
