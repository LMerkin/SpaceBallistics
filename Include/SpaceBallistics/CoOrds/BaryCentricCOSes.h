// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/BaryCentricCOSes.h":                 //
//                    The Inertial Co-Ordinate System                        //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"

namespace SpaceBallistics
{
  // Fwd Decls:
  class  TDB;

  //=========================================================================//
  // "BaryCEclCOS" Struct:                                                   //
  //=========================================================================//
  // BaryCentric Ecliptical COS:
  // Origin: Solar System BaryCenter
  // Axes  :
  //   X       : ICRS/BCRS ~ Dynamic Equinox  of J2000.0;
  //   XY plane: ICRS/BCSR ~ Dynamic Ecliptic of J2000.0;
  //             (ie NOT ICRS/BCRS itself, which is Equatorial, and NOT the La-
  //             place plane of the Solar System, which is the plane orthogonal
  //             to the angular momentum vector of the whole Solar System);
  //   Y, Z    : Derived from the above;
  // Obviously, using TDB as the associated TimeScale;
  // This definition is consistent with the conventions of the JPL Horizon eph-
  // emerides:
  //
  struct BaryCEclCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCEclCOS() = delete;  // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  template<Body B   = Body::UNDEFINED>
  using PosKVBarEcl = PosKV<BaryCEclCOS, B>;

  template<Body B   = Body::UNDEFINED>
  using VelKVBarEcl = VelKV<BaryCEclCOS, B>;

  // XXX: Currently no need to consider other Vectors in this COS yet...

  //=========================================================================//
  // "BaryCEqCOS" Struct:                                                    //
  //=========================================================================//
  // BaryCentric Equatorial COS:
  // As "BaryCEclCOS" above, but the XY plane is ~ the Mean Earth Equator of
  // J2000.0. Thus, this COS is the ICRS/BCRS itself:
  //
  struct BaryCEqCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCEqCOS() = delete;   // No objects construction at all!
  };

  // Alias:
  using BCRS = BaryCEqCOS;

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  template<Body B   = Body::UNDEFINED>
  using PosKVBarEq  = PosKV<BCRS, B>;

  template<Body B   = Body::UNDEFINED>
  using VelKVBarEq  = VelKV<BCRS, B>;

  // Aliases:
  template<Body B   = Body::UNDEFINED>
  using PosKV_BCRS  = PosKV<BCRS, B>;

  template<Body B   = Body::UNDEFINED>
  using VelKV_BCRS  = VelKV<BCRS, B>;

  // XXX: Currently no need to consider other Vectors in this COS yet...
}
// End namespace SpaceBallistics
