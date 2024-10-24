// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/BodyCentricCOSes.h":               //
//        Body-Centric Fixed (ICRF Axes) and Rotating Co-Ords Systems        //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "BodyCentricFixedCOS" Class:                                            //
  //=========================================================================//
  // Origin: Normally, the center of the Body Ellipsoid.
  // Axes  : X       : Equinox(2000.0);
  //         XY Plane: Body's Equator(2000.0);
  //         Z       : Body's North Pole
  // TimeScale       : TT for Earth, TDB otherwise
  // This, the axes of this COS are fixed in the ICRF, but the Origin is moving
  // in the inertial space, so this COS is not fully-inertial; it may sometimes
  // be considered "quasi-inertial".
  // NB    : This class stands for itself; no objects of it can be created
  //
  class TT;
  class TDB;

  template<Body BodyName>
  struct   BodyCentricFixedCOS
  {
    using  TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    BodyCentricFixedCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricFixedCOS"es:                  //
  //-------------------------------------------------------------------------//
  using GeoCentricFixedCOS     = BodyCentricFixedCOS<Body::Earth>;
  using SelenoCentricFixedCOS  = BodyCentricFixedCOS<Body::Moon>;
  using HelioCentricFixedCOS   = BodyCentricFixedCOS<Body::Sun>;
  using HermeoCentricFixedCOS  = BodyCentricFixedCOS<Body::Mercury>;
  using CytheroCentricFixedCOS = BodyCentricFixedCOS<Body::Venus>;
  using AreoCentricFixedCOS    = BodyCentricFixedCOS<Body::Mars>;
  using ZenoCentricFixedCOS    = BodyCentricFixedCOS<Body::Jupiter>;
  using CronoCentricFixedCOS   = BodyCentricFixedCOS<Body::Saturn>;
  using UranoCentricFixedCOS   = BodyCentricFixedCOS<Body::Uranus>;
  using PoseidoCentricFixedCOS = BodyCentricFixedCOS<Body::Neptune>;

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCentricFixedCOS":        //
  //-------------------------------------------------------------------------//
  // NB: Pos and Vel Vectors use "Len_km" ("K"), but other Vectors use "Len_m":
  //
  template<Body BodyName>
  using PosKVFix   = PosKV  <BodyCentricFixedCOS<BodyName>>;

  template<Body BodyName>
  using VelKVFix   = VelKV  <BodyCentricFixedCOS<BodyName>>;

  template<Body BodyName>
  using AccVFix    = AccV   <BodyCentricFixedCOS<BodyName>>;

  template<Body BodyName>
  using ForceVFix  = ForceV <BodyCentricFixedCOS<BodyName>>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...

  //=========================================================================//
  // "BodyCentricRotatingCOS" Class:                                         //
  //=========================================================================//
  // This COS is "embedded" in the Body and is rotating in the inertial space
  // along with the Body.
  // Origin: Normally, the center of the Body Ellipsoid.
  // Axes  : XY Plane: Body's current Equator;
  //         X       : To (lambda=0, phi=0) in the corresp BodyGraphic COS;
  //         Z       : Body's North Pole;
  // TimeScale       : Same convention as for BodyCentricFixedCOS.
  // This, the axes of this COS are rotating in the inertial (eg ICRF) space,
  // and the origin is moving as well, so it is definitely NON-INERTIAL.
  // in the inertial space, so this COS is not fully-inertial.
  // NB    : This class stands for itself; no objects of it cane be created :
  //
  template<Body BodyName>
  struct BodyCentricRotatingCOS
  {
    using  TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    BodyCentricRotatingCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Aliases: Proper Names of some "BodyCentricRotatingCOS"es:               //
  //-------------------------------------------------------------------------//
  using GeoCentricRotatingCOS     = BodyCentricRotatingCOS<Body::Earth>;
  using SelenoCentricRotatingCOS  = BodyCentricRotatingCOS<Body::Moon>;
  using HelioCentricRotatingCOS   = BodyCentricRotatingCOS<Body::Sun>;
  using HermeoCentricRotatingCOS  = BodyCentricRotatingCOS<Body::Mercury>;
  using CytheroCentricRotatingCOS = BodyCentricRotatingCOS<Body::Venus>;
  using AreoCentricRotatingCOS    = BodyCentricRotatingCOS<Body::Mars>;
  using ZenoCentricRotatingCOS    = BodyCentricRotatingCOS<Body::Jupiter>;
  using CronoCentricRotatingCOS   = BodyCentricRotatingCOS<Body::Saturn>;
  using UranoCentricRotatingCOS   = BodyCentricRotatingCOS<Body::Uranus>;
  using PoseidoCentricRotatingCOS = BodyCentricRotatingCOS<Body::Neptune>;

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in a "BodyCentricRotatingCOS":     //
  //-------------------------------------------------------------------------//
  // Again, using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BodyName>
  using PosKVRot   = PosKV  <BodyCentricRotatingCOS<BodyName>>;

  template<Body BodyName>
  using VelKVRot   = VelKV  <BodyCentricRotatingCOS<BodyName>>;

  template<Body BodyName>
  using AccVRot    = AccV   <BodyCentricRotatingCOS<BodyName>>;

  template<Body BodyName>
  using ForceVRot  = ForceV <BodyCentricRotatingCOS<BodyName>>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...
}
// End namespace SpaceBallistics
