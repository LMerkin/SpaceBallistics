// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/BodyCentricCOSes.h":               //
//           Body-Centric Fixed-Axes and Rotating Co-Ords Systems            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
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
  // Origin: Normally, the center of Body's Ellipsoid
  // Axes  : X       : Typically, Body's Vernal Equinox(J2000.0) = the Ascending
  //                   Node of the Sun orbit over Body's Equator  (ie Body's
  //                   Ecliptic) = the Descending Node of Body's  orbit over
  //                   Body's Equator
  //         XY Plane: Body's Equator(J2000.0)
  //         Z       : Body's North Pole
  // IN PARTICULAR, the axes of "GeoCentricEqFixCOS" are exactly those of ICRF!
  // TimeScale       : TT for Earth, TDB otherwise
  // This, the axes of this COS are fixed in the ICRF, but the Origin is moving
  // in the inertial space, so this COS is not fully-inertial
  // NB    : This class stands for itself; no objects of it can be created
  //
  template<Body BodyName>
  struct   BodyCentricEqFixCOS
  {
    using  TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    BodyCentricEqFixCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricEqFixCOS"es:                  //
  //-------------------------------------------------------------------------//
  using HelioCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Sun>;
  using HermeoCentricEqFixCOS  = BodyCentricEqFixCOS<Body::Mercury>;
  using CytheroCentricEqFixCOS = BodyCentricEqFixCOS<Body::Venus>;
  using GeoCentricEqFixCOS     = BodyCentricEqFixCOS<Body::Earth>;   // GCRS !!!
  using SelenoCentricEqFixCOS  = BodyCentricEqFixCOS<Body::Moon>;
  using AreoCentricEqFixCOS    = BodyCentricEqFixCOS<Body::Mars>;
  using ZenoCentricEqFixCOS    = BodyCentricEqFixCOS<Body::Jupiter>;
  using CronoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Saturn>;
  using UranoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::Uranus>;
  using PoseidoCentricEqFixCOS = BodyCentricEqFixCOS<Body::Neptune>;
  using HadeoCentricEqFixCOS   = BodyCentricEqFixCOS<Body::PlChB>;   // XXX ???

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCentricEqFixCOS":        //
  //-------------------------------------------------------------------------//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body BodyName>
  using PosKVEqFix      = PosKV <BodyCentricEqFixCOS<BodyName>>;

  template<Body BodyName>
  using VelKVEqFix      = VelKV <BodyCentricEqFixCOS<BodyName>>;

  template<Body BodyName>
  using AccVEqFix       = AccV  <BodyCentricEqFixCOS<BodyName>>;

  template<Body BodyName>
  using ForceVEqFix     = ForceV<BodyCentricEqFixCOS<BodyName>>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEqFixCOS":                        //
  //-------------------------------------------------------------------------//
  using PosKVGeoEqFix   = PosKVEqFix <Body::Earth>;
  using VelKVGeoEqFix   = VelKVEqFix <Body::Earth>;
  using AccVGeoEqFix    = AccVEqFix  <Body::Earth>;
  using ForceVGeoEqFix  = ForceVEqFix<Body::Earth>;

  //=========================================================================//
  // "BodyCentricEclFixCOS":                                                 //
  //=========================================================================//
  // BodyCentric Ecliptical Fixed-Axes COS:
  // Origin: Normally, the center of Body's Ellipsoid
  // Axes  : X       : Same as the X axis of the corresp "BodyCentricEqFixCOS"
  //         XY Plane: Body's Mean Orbital Plane(J2000.0)
  //         Z       : "North Orbital Pole"
  // TimeScale       : TT for Earth, TDB otherwise
  // Again, the axes of this COS are fixed in the ICRF, but the Origin is moving
  // in the inertial space, so this COS is also not fully-inertial
  // NB    : This class stands for itself; no objects of it can be created
  //
  // Similar to "BodyCentricEqFixCOS" above, but the XY plane is Body's Orbital
  // Plane (rather than Body's Equator) @ J2000.0. XXX: Again, is this convent-
  // ion optimal for extraterrestrial Bodies?
  //
  template<Body BodyName>
  struct   BodyCentricEclFixCOS
  {
    using  TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

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
  template<Body BodyName>
  using PosKVEclFix     = PosKV <BodyCentricEclFixCOS<BodyName>>;

  template<Body BodyName>
  using VelKVEclFix     = VelKV <BodyCentricEclFixCOS<BodyName>>;

  template<Body BodyName>
  using AccVEclFix      = AccV  <BodyCentricEclFixCOS<BodyName>>;

  template<Body BodyName>
  using ForceVEclFix    = ForceV<BodyCentricEclFixCOS<BodyName>>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEclFixCOS":                       //
  //-------------------------------------------------------------------------//
  using PosKVGeoEclFix  = PosKVEclFix <Body::Earth>;
  using VelKVGeoEclFix  = VelKVEclFix <Body::Earth>;
  using AccVGeoEclFix   = AccVEclFix  <Body::Earth>;
  using ForceVGeoEclFix = ForceVEclFix<Body::Earth>;

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "HelioCentricEclFixCOS":                     //
  //-------------------------------------------------------------------------//
  using PosKVHelEclFix  = PosKVEclFix <Body::Sun>;
  using VelKVHelEclFix  = VelKVEclFix <Body::Sun>;
  using AccVHelEclFix   = AccVEclFix  <Body::Sun>;
  using ForceVHelEclFix = ForceVEclFix<Body::Sun>;

  //=========================================================================//
  // "BodyCentricRotCOS" Class:                                              //
  //=========================================================================//
  // BodyCentric Rotating-Axes COS (obviously Equatorial):
  // This COS is "embedded" in the Body and is rotating in the inertial space
  // along with the Body.
  // Origin: Normally, the center of the Body Ellipsoid.
  // Axes  : XY Plane: Body's CURRENT Equator (NOT the J2000.0 Equator!)
  //         X       : To (lambda=0, phi=0) in the corresp BodyGraphic COS
  //         Z       : Body's North Pole
  // TimeScale       : Same convention as for BodyCentric{Eq,Ecl}FixCOS.
  // This, the axes of this COS are rotating in the inertial (eg ICRF) space,
  // the rotation axis may precess  in the inertial space,  and the origin is
  // moving as well, so this system is STRONGLY NON-INERTIAL.
  // NB    : This class stands for itself; no objects of it cane be created :
  //
  template<Body BodyName>
  struct BodyCentricRotCOS
  {
    using  TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    BodyCentricRotCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricRotCOS"es:                    //
  //-------------------------------------------------------------------------//
  using HelioCentricRotCOS   = BodyCentricRotCOS<Body::Sun>;
  using HermeoCentricRotCOS  = BodyCentricRotCOS<Body::Mercury>;
  using CytheroCentricRotCOS = BodyCentricRotCOS<Body::Venus>;
  using GeoCentricRotCOS     = BodyCentricRotCOS<Body::Earth>;    // ITRS !!!
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
  template<Body BodyName>
  using PosKVRot     = PosKV <BodyCentricRotCOS<BodyName>>;

  template<Body BodyName>
  using VelKVRot     = VelKV <BodyCentricRotCOS<BodyName>>;

  template<Body BodyName>
  using AccVRot      = AccV  <BodyCentricRotCOS<BodyName>>;

  template<Body BodyName>
  using ForceVRot    = ForceV<BodyCentricRotCOS<BodyName>>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCentricEqRotCOS":                        //
  //-------------------------------------------------------------------------//
  using PosKVGeoRot  = PosKVRot <Body::Earth>;
  using VelKVGeoRot  = VelKVRot <Body::Earth>;
  using AccVGeoRot   = AccVRot  <Body::Earth>;
  using ForceVGeoRot = ForceVRot<Body::Earth>;
}
// End namespace SpaceBallistics
