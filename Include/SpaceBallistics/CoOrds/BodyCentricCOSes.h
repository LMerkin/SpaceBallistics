// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/BodyCentricCOSes.h":               //
//           Body-Centric Fixed-Axes and Rotating Co-Ords Systems            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include <type_traits>

namespace SpaceBallistics
{
  //-------------------------------------------------------------------------//
  // Fwd Decls of Time Scales:                                               //
  //-------------------------------------------------------------------------//
  class TT;
  class TDB;

  //=========================================================================//
  // "BodyCentricEqFixCOS" Class:                                            //
  //=========================================================================//
  // BodyCentric Equatorial Fixed-Axes COS:
  // Origin: Normally, the center of Body's Ellipsoid;
  // Axes  : X       : Typically, ~ Body's Dynamic Vernal Equinox(J2000.0)  =
  //                   the Ascending Node of the Sun orbit over Body's Equator
  //                   (ie "Body's Ecliptic") = the Descending Node of Body's
  //                   orbit over Body's Equator;
  //         XY Plane: ~ Body's Mean Equator(J2000.0);
  //         Z       : ~ Body's North Pole;
  // TimeScale       : TT for Earth, TDB otherwise (FIXME: It would be better
  //                   to provide Body-centric relativistic TimeScales for all
  //                   Bodies, not just for Earth).
  // For Earth, this is exactly the GCRS!
  // This, the axes of this COS are fixed in the ICRS, but the Origin is moving
  // in the inertial space, so this COS is not fully-inertial:
  //
  template<Body Origin>
  struct   BodyCentricEqFixCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = true;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // This struct stands for itself; no objects of it can be created:
    BodyCentricEqFixCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricEqFixCOS"es:                  //
  //-------------------------------------------------------------------------//
  using HelioCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Sun>;
  using HermeoCentricEqFixCOS  = BodyCentricEqFixCOS<Body::Mercury>;
  using CytheroCentricEqFixCOS = BodyCentricEqFixCOS<Body::Venus>;
  using GeoCentricEqFixCOS     = BodyCentricEqFixCOS<Body::Earth>;
  using GCRS                   = GeoCentricEqFixCOS;
  using SelenoCentricEqFixCOS  = BodyCentricEqFixCOS<Body::Moon>;
  using AreoCentricEqFixCOS    = BodyCentricEqFixCOS<Body::Mars>;
  using ZenoCentricEqFixCOS    = BodyCentricEqFixCOS<Body::Jupiter>;
  using CronoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Saturn>;
  using UranoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Uranus>;
  using PoseidoCentricEqFixCOS = BodyCentricEqFixCOS<Body::Neptune>;
  using HadeoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::PlChB>;   // XXX ???

  //=========================================================================//
  // Position, Velocity and other Vectors in a "BodyCentricEqFixCOS":        //
  //=========================================================================//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVEqFix  = PosKV <BodyCentricEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVEqFix  = VelKV <BodyCentricEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVEqFix   = AccV  <BodyCentricEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVEqFix = ForceV<BodyCentricEqFixCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEqFixCOS":                        //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using PosKV_GCRS  = PosKVEqFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKV_GCRS  = VelKVEqFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccV_GCRS   = AccVEqFix  <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceV_GCRS = ForceVEqFix<Body::Earth, B>;

  //-------------------------------------------------------------------------//
  // Translation of Origins between ICRS/BCRS and GCRS:                      //
  //-------------------------------------------------------------------------//
  // The orientation of axes is the same, so:
  // (Bary,  Earth) + (Earth, Body) => (Bary, Body ):
  //
  template<typename DQ,     Body B>
  Vector3D<DQ, BaryCentricEqCOS, B> operator+
  (
    Vector3D<DQ, BaryCentricEqCOS,   Body::Earth> const& a_earth,
    Vector3D<DQ, GeoCentricEqFixCOS, B>           const& a_geo
  )
  {
    return Vector3D<DQ, BaryCentricEqCOS, B>
           { a_earth.x() + a_geo.x(),
             a_earth.y() + a_geo.y(),
             a_earth.z() + a_geo.z()
           };
  }

  // (Bary,  Body) - (Bary, Earth) = (Earth, Body):
  //
  template<typename DQ,      Body  B>
  Vector3D<DQ, GeoCentricEqFixCOS, B> operator-
  (
    Vector3D<DQ, BaryCentricEqCOS, B>           const& a_body,
    Vector3D<DQ, BaryCentricEqCOS, Body::Earth> const& a_earth
  )
  {
    return Vector3D<DQ, GeoCentricEqFixCOS, B>
           { a_body.x() - a_earth.x(),
             a_body.y() - a_earth.y(),
             a_body.z() - a_earth.z()
           };
  }

  //=========================================================================//
  // "BodyCentricEclFixCOS":                                                 //
  //=========================================================================//
  // BodyCentric Ecliptical Fixed-Axes COS:
  // Origin: Normally, the center of Body's Ellipsoid;
  // Axes  : X       : Same as the X axis of the corresp "BodyCentricEqFixCOS"
  //         XY Plane: Body's Mean Orbital Plane(J2000.0); for Earth, Ecliptic
  //                   is defined as the Mean Earth Orbital Plane;
  //         Z       : "North Orbital Pole";
  // TimeScale       : TT for Earth, TDB otherwise (FIXME: again, it would be
  //                   better to provide relativistic Body-centric TimeScales
  //                   for all Bodies, not just for Earth);
  // Again, the axes of this COS are fixed in the ICRS, but the Origin is moving
  // in the inertial space, so this COS is also not fully-inertial
  //
  // Similar to "BodyCentricEqFixCOS" above, but the XY plane is ~ Body's Dyna-
  // mic Orbital Plane (rather than Body's Equator) @ J2000.0:
  //
  template<Body Origin>
  struct   BodyCentricEclFixCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = true;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // This struct stands for itself; no objects of it can be created:
    BodyCentricEclFixCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricEclFixCOS"es:                 //
  //-------------------------------------------------------------------------//
  using HelioCentricEclFixCOS   = BodyCentricEclFixCOS<Body::Sun>;
  using HermeoCentricEclFixCOS  = BodyCentricEclFixCOS<Body::Mercury>;
  using CytheroCentricEclFixCOS = BodyCentricEclFixCOS<Body::Venus>;
  using GeoCentricEclFixCOS     = BodyCentricEclFixCOS<Body::Earth>;
  using SelenoCentricEclFixCOS  = BodyCentricEclFixCOS<Body::Moon>;
  using AreoCentricEclFixCOS    = BodyCentricEclFixCOS<Body::Mars>;
  using ZenoCentricEclFixCOS    = BodyCentricEclFixCOS<Body::Jupiter>;
  using CronoCentricEclFixCOS   = BodyCentricEclFixCOS<Body::Saturn>;
  using UranoCentricEclFixCOS   = BodyCentricEclFixCOS<Body::Uranus>;
  using PoseidoCentricEclFixCOS = BodyCentricEclFixCOS<Body::Neptune>;
  using HadeoCentricEclFixCOS   = BodyCentricEclFixCOS<Body::PlChB>; // XXX ???

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCentricEclFixCOS":       //
  //-------------------------------------------------------------------------//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVEclFix     = PosKV <BodyCentricEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVEclFix     = VelKV <BodyCentricEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVEclFix      = AccV  <BodyCentricEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVEclFix    = ForceV<BodyCentricEclFixCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEclFixCOS":                       //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using PosKVGeoEclFix  = PosKVEclFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKVGeoEclFix  = VelKVEclFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccVGeoEclFix   = AccVEclFix  <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceVGeoEclFix = ForceVEclFix<Body::Earth, B>;

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "HelioCentricEclFixCOS":                     //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using PosKVHelEclFix  = PosKVEclFix <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using VelKVHelEclFix  = VelKVEclFix <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using AccVHelEclFix   = AccVEclFix  <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using ForceVHelEclFix = ForceVEclFix<Body::Sun, B>;

  //-------------------------------------------------------------------------//
  // Translation of Origins between Ecliptical {Bari,Geo}Centric COSes:      //
  //-------------------------------------------------------------------------//
  // The orientation of axes is the same, so:
  // (Bary,  Earth) + (Earth, Body) => (Bary, Body ):
  //
  template<typename DQ,      Body B>
  Vector3D<DQ, BaryCentricEclCOS, B> operator+
  (
    Vector3D<DQ, BaryCentricEclCOS,   Body::Earth> const& a_earth,
    Vector3D<DQ, GeoCentricEclFixCOS, B>           const& a_geo
  )
  {
    return Vector3D<DQ, BaryCentricEclCOS, B>
           { a_earth.x() + a_geo.x(),
             a_earth.y() + a_geo.y(),
             a_earth.z() + a_geo.z()
           };
  }

  // (Bary,  Body) - (Bary, Earth) = (Earth, Body):
  //
  template<typename DQ,       Body  B>
  Vector3D<DQ, GeoCentricEclFixCOS, B> operator-
  (
    Vector3D<DQ, BaryCentricEclCOS, B>           const& a_body,
    Vector3D<DQ, BaryCentricEclCOS, Body::Earth> const& a_earth
  )
  {
    return Vector3D<DQ, GeoCentricEclFixCOS, B>
           { a_body.x() - a_earth.x(),
             a_body.y() - a_earth.y(),
             a_body.z() - a_earth.z()
           };
  }

  //=========================================================================//
  // "BodyCentricRotCOS" Class:                                              //
  //=========================================================================//
  // BodyCentric Rotating-Axes COS (obviously Equatorial):
  // This COS is "embedded" in the Body and is rotating in the inertial space
  // along with the Body.
  // Origin: Normally, the center of the Body Ellipsoid;
  // Axes  : XY Plane: Body's CURRENT Dynamical Equator (NOT the J2000.0 Equator
  //                   and NOT the Mean Equator of the date!);
  //         X       : To (lambda=0, phi=0) in the corresp BodyGraphic COS;
  //         Z       : Body's North Pole;
  // TimeScale       : Same convention as for BodyCentric{Eq,Ecl}FixCOS.
  // This, the axes of this COS are rotating in the inertial (eg ICRS) space,
  // the rotation axis may precess  in the inertial space,  and the origin is
  // moving as well, so this system is STRONGLY NON-INERTIAL.
  // For Earth, this is exactly the ITRS!
  //
  template<Body Origin>
  struct BodyCentricRotCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = false;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // This struct stands for itself; no objects of it can be created:
    BodyCentricRotCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricRotCOS"es:                    //
  //-------------------------------------------------------------------------//
  using HelioCentricRotCOS   = BodyCentricRotCOS<Body::Sun>;
  using HermeoCentricRotCOS  = BodyCentricRotCOS<Body::Mercury>;
  using CytheroCentricRotCOS = BodyCentricRotCOS<Body::Venus>;
  using GeoCentricRotCOS     = BodyCentricRotCOS<Body::Earth>;
  using ITRS                 = GeoCentricRotCOS;
  using SelenoCentricRotCOS  = BodyCentricRotCOS<Body::Moon>;
  using AreoCentricRotCOS    = BodyCentricRotCOS<Body::Mars>;
  using ZenoCentricRotCOS    = BodyCentricRotCOS<Body::Jupiter>;
  using CronoCentricRotCOS   = BodyCentricRotCOS<Body::Saturn>;
  using UranoCentricRotCOS   = BodyCentricRotCOS<Body::Uranus>;
  using PoseidoCentricRotCOS = BodyCentricRotCOS<Body::Neptune>;
  using HadeoCentricRotCOS   = BodyCentricRotCOS<Body::PlChB>;    // XXX ???

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCentricRotCOS":          //
  //-------------------------------------------------------------------------//
  // Again, using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVRot    = PosKV <BodyCentricRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVRot    = VelKV <BodyCentricRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVRot     = AccV  <BodyCentricRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVRot   = ForceV<BodyCentricRotCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEqRotCOS":                        //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using PosKV_ITRS  = PosKVRot <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKV_ITRS  = VelKVRot <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccV_ITRS   = AccVRot  <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceV_ITRS = ForceVRot<Body::Earth, B>;
}
// End namespace SpaceBallistics
