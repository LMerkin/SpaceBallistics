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
  // "BaryCentricEclCOS" Struct:                                             //
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
  struct BaryCentricEclCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCentricEclCOS() = delete;  // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  template<Body B   = Body::UNDEFINED>
  using PosKVBarEcl = PosKV<BaryCentricEclCOS, B>;

  template<Body B   = Body::UNDEFINED>
  using VelKVBarEcl = VelKV<BaryCentricEclCOS, B>;

  // XXX: Currently no need to consider other Vectors in this COS yet...

  //=========================================================================//
  // "BaryCentricEqCOS" Struct:                                              //
  //=========================================================================//
  // BaryCentric Equatorial COS:
  // As "BaryCentricEclCOS" above, but the XY plane is ~ the Mean Earth Equator
  // of J2000.0. Thus, this COS is the ICRS/BCRS itself:
  //
  struct BaryCentricEqCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCentricEqCOS() = delete;   // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  template<Body B   = Body::UNDEFINED>
  using PosKVBarEq  = PosKV<BaryCentricEqCOS, B>;

  template<Body B   = Body::UNDEFINED>
  using VelKVBarEq  = VelKV<BaryCentricEqCOS, B>;

  // XXX: Currently no need to consider other Vectors in this COS yet...
}
// End namespace SpaceBallistics
