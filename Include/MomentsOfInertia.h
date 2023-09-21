// vim:ts=2:et
//===========================================================================//
//                            "MomentsOfInertia.h":                          //
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

  //=========================================================================//
  // "MoIS" Struct:                                                          //
  //=========================================================================//
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
  // MoI Computation for Thin Shells:                                        //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "MoIS_TrCone_XY":                                                       //
  //-------------------------------------------------------------------------//
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
  // Result  : {Moment of Inertia wrt the OY axis per Surface Density [L^4],
  //            Side Surface Area [L^2]}:
  //
  MoIS MoIS_TrCone_XY
  (
    Len_m   a_d0,     // Base diameter at a_x0
    Len_m   a_d1,     // Base diameter at the other section (at x = x1 >= x0)
    Len_m   a_h,      // Must be > 0
    Len_m   a_x0,
    double  a_alpha   // Dimension-less, in [-Pi/2 .. Pi/2]
  );

  //-------------------------------------------------------------------------//
  // "MoIS_TrCone_XZ":                                                       //
  //-------------------------------------------------------------------------//
  // The geometry is similar to that in "MoIS_TrConeXY", but the axis of the
  // truncated cone is now located in the OXZ plane, at the angle "alpha" to
  // the OX axis (|alpha| <= Pi/2);  (x0,z0) are the coords of the "initial"
  // point of the axis segment (with the lowest X-coord), corresp to the base
  // diameter "d0"; "d1" is the other base's diameter. NB: in this case, both
  // "x0" and "z0" are required. The result is invariant under the transform
  // z0 <-> (-z0), alpha <-> (-alpha):
  //
  MoIS MoIS_TrCone_XZ
  (
    Len_m   a_d0,     // Base diameter at (a_x0, a_z0)
    Len_m   a_d1,     // Base diameter at the other section (at x = x1 >= x0)
    Len_m   a_h,      // Must be > 0
    Len_m   a_x0,
    Len_m   a_z0,
    double  a_alpha   // Dimension-less, in [-Pi/2 .. Pi/2]
  );

  //-------------------------------------------------------------------------//
  // "MoIS_SpherSegm_XY":                                                    //
  //-------------------------------------------------------------------------//
  // Similar to "MoIS_TrCone_XY", but for a Spherical Segment of the base diam-
  // eter "d" and the height (from the base plane to the pole) "h", where we
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
  // "is_pos_facing" flag:
  //
  MoIS MoIS_SpherSegm_XY
  (
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    bool    a_is_pos_facing
  );

  //-------------------------------------------------------------------------//
  // "MoIS_SpherSegm_XZ":                                                    //
  //-------------------------------------------------------------------------//
  // Similar to "MoIS_SpherSegm_XY", but the axis of the spherical segment is
  // in the OXZ plane, at the angle "alpha" to the positive direction  of the
  // OX axis, and the coords of the base center in that plane are (x0, z0):
  //
  MoIS MoIS_SpherSegm_XZ
  (
    Len_m   a_d,            // Base diameter
    Len_m   a_h,            // Height: 0 < a_h <= a_d/2.0
    Len_m   a_x0,           // X-coord of the base center (NOT of the pole!)
    Len_m   a_z0,           // Z-coord of the base center (NOT of the pole!)
    double  a_alpha,        // Dimension-less, in [-Pi/2 .. Pi/2]
    bool    a_is_pos_facing
  );
}
// End namespace SpaceBallistics
