// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/BodyCentricCOSes.h":               //
//           Body-Centric Fixed-Axes and Rotating Co-Ords Systems            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // IMPORTANT REMARK:                                                       //
  //=========================================================================//
  // For the avoidance of doubt, all BodyCentric COSes are NON-INERTIAL.  They
  // are considered to be SnapShots; in particular, the velocity, acceleration
  // etc of the Body itself in any BodyCentric COS are NOT identical 0s.
  // The transformations between BodyCentric COSes and other COSes (eg the in-
  // ertial BaryCentric ones) would necessarily involve the SnapShot time as a
  // parameter. Ideally, those parameters would be installed as properties in
  // the BodyCentric COS types; unfortunately,  this is not possible because
  // the SnapShot TimeStamp is usually known at run-time only.   Rather, the
  // SnapShot TimeStamp of the SnapShot is installed in all "Vector3D"s refer-
  // ring to any COS.
  //
  //=========================================================================//
  // "BodyCEqFixCOS" Class:                                                  //
  //=========================================================================//
  // BodyCentric Equatorial Fixed-Axes COS:
  // Origin: Normally, the center of Body's Ellipsoid;
  // Axes  : X       : Typically, ~ Body's Dynamic Vernal Equinox(J2000.0)  =
  //                   the Ascending Node of the Sun orbit over Body's Equator
  //                   (ie "Body's Ecliptic") = the Descending Node of Body's
  //                   orbit over Body's Equator;
  //         XY Plane: ~ Body's Mean Equator(J2000.0);
  //         Z       : ~ Body's Mean North Pole;
  // TimeScale       : TT for Earth, TDB otherwise (XXX: It would be better to
  //                   provide Body-centric relativistic TimeScales for all Bo-
  //                   dies, not just for Earth).
  // For Earth, this is exactly the GCRS! That is, the Earth Mean Equator (w/o
  // short-term Nutations) and the Dynamic Ecliptic and thus Equinox of J2000.0
  // are used.
  // This, the axes of this COS are fixed in the ICRS, but the Origin is moving
  // in the inertial space, so this COS is not fully-inertial:
  //
  template<Body Origin>
  struct   BodyCEqFixCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = true;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // NB: The Origin MUST NOT be UnDefined:
    static_assert(Origin != Body::UNDEFINED);

    // This struct stands for itself; no objects of it can be created:
    BodyCEqFixCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCEqFixCOS"es:                        //
  //-------------------------------------------------------------------------//
  using HelioCEqFixCOS   = BodyCEqFixCOS<Body::Sun>;
  using HermeoCEqFixCOS  = BodyCEqFixCOS<Body::Mercury>;
  using CytheroCEqFixCOS = BodyCEqFixCOS<Body::Venus>;
  using GeoCEqFixCOS     = BodyCEqFixCOS<Body::Earth>;
  using GCRS             = GeoCEqFixCOS;
  using SelenoCEqFixCOS  = BodyCEqFixCOS<Body::Moon>;
  using AreoCEqFixCOS    = BodyCEqFixCOS<Body::Mars>;
  using ZenoCEqFixCOS    = BodyCEqFixCOS<Body::Jupiter>;
  using CronoCEqFixCOS   = BodyCEqFixCOS<Body::Saturn>;
  using UranoCEqFixCOS   = BodyCEqFixCOS<Body::Uranus>;
  using PoseidoCEqFixCOS = BodyCEqFixCOS<Body::Neptune>;
  using HadeoCEqFixCOS   = BodyCEqFixCOS<Body::PlChB>;   // XXX ???

  //=========================================================================//
  // Position, Velocity and other Vectors in a "BodyCEqFixCOS":              //
  //=========================================================================//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using DimLessVEqFix = DimLessV<BodyCEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVEqFix    = PosKV   <BodyCEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVEqFix    = VelKV   <BodyCEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVEqFix     = AccV    <BodyCEqFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVEqFix   = ForceV  <BodyCEqFixCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCEqFixCOS":                              //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using DimLessV_GCRS = DimLessVEqFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using PosKV_GCRS    = PosKVEqFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKV_GCRS    = VelKVEqFix <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccV_GCRS     = AccVEqFix  <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceV_GCRS   = ForceVEqFix<Body::Earth, B>;

  //=========================================================================//
  // "BodyCEclFixCOS":                                                       //
  //=========================================================================//
  // BodyCentric Ecliptical Fixed-Axes COS:
  // Origin: Normally, the center of Body's Ellipsoid;
  // Axes  : X       : Same as the X axis of the corresp "BodyCEqFixCOS"
  //         XY Plane: Body's Dynamic Orbital Plane(J2000.0); for Earth, it is
  //                   the Dynamic Ecliptic of J2000.0;
  //         Z       : Dynamic "North Orbital Pole";
  // TimeScale       : TT for Earth, TDB otherwise (XXX: again, it would be bet-
  //                   ter to provide relativistic Body-centric TimeScales  for
  //                   all Bodies, not just for Earth);
  // Again, the axes of this COS are fixed in the ICRS, but the Origin is moving
  // in the inertial space, so this COS is also not fully-inertial
  //
  // Similar to "BodyCEqFixCOS" above, but the XY plane is Body's Dynamic Orbit-
  // al Plane (rather than the Mean Equator) @ J2000.0. Ie, for Earth, the XY
  // Plane is the Dynamic Ecliptic of J2000.0:
  //
  template<Body Origin>
  struct   BodyCEclFixCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = true;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // NB: The Origin MUST NOT be UnDefined:
    static_assert(Origin != Body::UNDEFINED);

    // This struct stands for itself; no objects of it can be created:
    BodyCEclFixCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCEclFixCOS"es:                       //
  //-------------------------------------------------------------------------//
  using HelioCEclFixCOS   = BodyCEclFixCOS<Body::Sun>;
  using HermeoCEclFixCOS  = BodyCEclFixCOS<Body::Mercury>;
  using CytheroCEclFixCOS = BodyCEclFixCOS<Body::Venus>;
  using GeoCEclFixCOS     = BodyCEclFixCOS<Body::Earth>;
  using SelenoCEclFixCOS  = BodyCEclFixCOS<Body::Moon>;
  using AreoCEclFixCOS    = BodyCEclFixCOS<Body::Mars>;
  using ZenoCEclFixCOS    = BodyCEclFixCOS<Body::Jupiter>;
  using CronoCEclFixCOS   = BodyCEclFixCOS<Body::Saturn>;
  using UranoCEclFixCOS   = BodyCEclFixCOS<Body::Uranus>;
  using PoseidoCEclFixCOS = BodyCEclFixCOS<Body::Neptune>;
  using HadeoCEclFixCOS   = BodyCEclFixCOS<Body::PlChB>; // XXX ???

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCEclFixCOS":             //
  //-------------------------------------------------------------------------//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using DimLessVEclFix  = DimLessV<BodyCEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVEclFix     = PosKV   <BodyCEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVEclFix     = VelKV   <BodyCEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVEclFix      = AccV    <BodyCEclFixCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVEclFix    = ForceV  <BodyCEclFixCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCEclFixCOS":                             //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using DimLessVGeoEclFix = DimLessVEclFix<Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using PosKVGeoEclFix    = PosKVEclFix   <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKVGeoEclFix    = VelKVEclFix   <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccVGeoEclFix     = AccVEclFix    <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceVGeoEclFix   = ForceVEclFix  <Body::Earth, B>;

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "HelioCEclFixCOS":                           //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using DimLessVHelEclFix = DimLessVEclFix<Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using PosKVHelEclFix    = PosKVEclFix <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using VelKVHelEclFix    = VelKVEclFix <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using AccVHelEclFix     = AccVEclFix  <Body::Sun, B>;

  template<Body B = Body::UNDEFINED>
  using ForceVHelEclFix   = ForceVEclFix<Body::Sun, B>;

  //=========================================================================//
  // "BodyCRotCOS" Class:                                                    //
  //=========================================================================//
  // BodyCentric Rotating-Axes COS (obviously Equatorial):
  // This COS is "embedded" in the Body and is rotating in the inertial space
  // along with the Body.
  // Origin: Normally, the center of the Body Ellipsoid;
  // Axes  : XY Plane: Body's CURRENT Dynamical Equator (NOT the J2000.0  Mean
  //                   Equator and NOT the Mean Equator of the date!); NB: the
  //                   difference between the Mean and the Dynamic Equator is
  //                   in accounting of short-term Nutations (no / yes, resp);
  //         X       : To (lambda=0, phi=0) in the corresp BodyGraphic COS;
  //         Z       : Body's Dynamic North Pole;
  // TimeScale       : Same convention as for "BodyC{Eq,Ecl}FixCOS".
  // This, the axes of this COS are rotating in the inertial (eg ICRS) space,
  // the rotation axis may precess  in the inertial space,  and the origin is
  // moving as well, so this system is STRONGLY NON-INERTIAL.
  // For Earth, this is exactly the ITRS!
  //
  template<Body Origin>
  struct BodyCRotCOS
  {
    constexpr static   Body BaseBody       = Origin;
    constexpr static   bool HasFixedAxes   = false;
    constexpr static   bool HasFixedOrigin = false;
    using  TimeScale = std::conditional_t<Origin == Body::Earth, TT, TDB>;

    // NB: The Origin MUST NOT be UnDefined:
    static_assert(Origin != Body::UNDEFINED);

    // This struct stands for itself; no objects of it can be created:
    BodyCRotCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCRotCOS"es:                          //
  //-------------------------------------------------------------------------//
  using HelioCRotCOS   = BodyCRotCOS<Body::Sun>;
  using HermeoCRotCOS  = BodyCRotCOS<Body::Mercury>;
  using CytheroCRotCOS = BodyCRotCOS<Body::Venus>;
  using GeoCRotCOS     = BodyCRotCOS<Body::Earth>;
  using ITRS           = GeoCRotCOS;
  using SelenoCRotCOS  = BodyCRotCOS<Body::Moon>;
  using AreoCRotCOS    = BodyCRotCOS<Body::Mars>;
  using ZenoCRotCOS    = BodyCRotCOS<Body::Jupiter>;
  using CronoCRotCOS   = BodyCRotCOS<Body::Saturn>;
  using UranoCRotCOS   = BodyCRotCOS<Body::Uranus>;
  using PoseidoCRotCOS = BodyCRotCOS<Body::Neptune>;
  using HadeoCRotCOS   = BodyCRotCOS<Body::PlChB>;    // XXX ???

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCRotCOS":                //
  //-------------------------------------------------------------------------//
  // Again, using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body Origin, Body B = Body::UNDEFINED>
  using DimLessVRot   = DimLessV<BodyCRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using PosKVRot      = PosKV   <BodyCRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using VelKVRot      = VelKV   <BodyCRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using AccVRot       = AccV    <BodyCRotCOS<Origin>, B>;

  template<Body Origin, Body B = Body::UNDEFINED>
  using ForceVRot     = ForceV  <BodyCRotCOS<Origin>, B>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //-------------------------------------------------------------------------//
  // Aliases for Vectors in the "GeoCEqRotCOS":                              //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using DimLessV_ITRS = DimLessVRot<Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using PosKV_ITRS    = PosKVRot   <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using VelKV_ITRS    = VelKVRot   <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using AccV_ITRS     = AccVRot    <Body::Earth, B>;

  template<Body B = Body::UNDEFINED>
  using ForceV_ITRS   = ForceVRot  <Body::Earth, B>;
}
// End namespace SpaceBallistics
