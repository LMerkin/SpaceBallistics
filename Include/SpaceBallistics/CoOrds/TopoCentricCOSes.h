// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/TopoCentrticCOSes.h":               //
//                      Body-TopoCentric Co-Ord Systems                      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/Locations.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // Fwd Decls of TimeScales:                                                //
  //=========================================================================//
  class TT;
  class TDB;

  //=========================================================================//
  // "TopoCentricEqFixCOS" Struct:                                           //
  //=========================================================================//
  // The Axes are as in "BodyCentricEqFixCOS" (in particular, for Earth they are
  // ICRS/GCRS), but the Origin is a Location ("L") on the Body,  NOT the Body
  // Center.
  // Suitable for computation of TopoCentric Ephemerides in astronomical appls.
  // The Origin of this COS is rotating (and otherwise moving) in the inertial
  // space, but the Axes are fixed:
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCentricEqFixCOS
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = true;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCentricEqFixCOS() = delete; // No objs of this struct are to be created
  };

  // In particular, for Earth we call it "TopCentricGCRS":
  template<Location<Body::Earth> const* L>
  using TopoCentricGCRS = TopoCentricEqFixCOS<Body::Earth, L>;

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors in this COS:                              //
  //-------------------------------------------------------------------------//
  // We probably do not need any other DQs here:
  //
  template<Body BBody,  Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopEqFix = PosKV<TopoCentricEqFixCOS<BBody, L>, B>;

  template<Body BBody,  Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopEqFix = VelKV<TopoCentricEqFixCOS<BBody, L>, B>;

  // Aliases for Vectors in the "TopoCentricGCRS":
  template<Location<Body::Earth> const* L,        Body B = Body::UNDEFINED>
  using PosKVTopGCRS  = PosKV<TopoCentricGCRS<L>, B>;

  template<Location<Body::Earth> const* L,        Body B = Body::UNDEFINED>
  using VelKVTopGCRS  = VelKV<TopoCentricGCRS<L>, B>;

  //=========================================================================//
  // "TopoCentricRotCOS_ENZ" Struct:                                         //
  //=========================================================================//
  // Origin   : A point on the Body Surface given by "L";
  // Axes     : X=LocalEast   (tangential to the Parallel);
  //            Y=LocalNorth  (tangential to the Meridian);
  //            Z=LocalZenith (normal to the Ellipsoid Surface);
  // TimeScale: TT for Earth, TDB otherwise;
  // This COS is HIGHLY NON-INERTIAL: its Origin is rotating (and otherwise mov-
  // ing) in the inertial space, and the axes are rotating as well.
  // Suitable eg for modeling of Launch Ops and for Ground Tracks!
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCentricRotCOS_ENZ
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCentricRotCOS_ENZ() = delete;
    // No objs of this struct are to be created
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others;
  // the "ENZ" suffix stands for "East(X)-North(Y)-Zenith(Z)":
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopRotENZ     = PosKV   <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopRotENZ     = VelKV   <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopRotENZ      = VelV    <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopRotENZ      = AccV    <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopRotENZ    = ForceV  <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopRotENZ   = AngVelV <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopRot_ENZ  = AngAccV <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopRot_ENZ  = AngMomV <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopRot_ENZ    = TorqV   <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopRot_ENZ     = MoIV    <TopoCentricRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopRot_ENZ = MoIRateV<TopoCentricRotCOS_ENZ<BBody, L>, B>;

  //=========================================================================//
  // "TopoCentricRotCOS" Struct:                                             //
  //=========================================================================//
  // Origin   : A point on the Body Surface given by "L";
  // XYZ Axes : Parallel to those of the "BodyCentricRotCOS"
  //            (ie for Earth they are parallel to the ITRS axes);
  // TimeScale: TT for Earth, TDB otherwise;
  // This COS is HIGHLY NON-INERTIAL: its Origin is rotating (and otherwise mov-
  // ing) in the inertial space, and the axes are rotating as well.
  // Suitable eg for modeling of Launch Ops and for Ground Tracks!
  // XXX: This COS (and its Vectors) do not have any suffix, unlike the above
  // "TopoCentricRotCOS_ENZ":
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCentricRotCOS
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCentricRotCOS() = delete;
    // No objs of this struct are to be created
  };

  //-----------------------------------------------------------------------//
  // Technical Stuff, for Template Meta-Programming:                       //
  //-----------------------------------------------------------------------//
  // "IsTopoCentricRotCOS":
  template<typename COS>
  constexpr inline bool IsTopoCentricRotCOS = false;

  template<Body BBody,  Location<BBody> const* L>
  constexpr inline bool IsTopoCentricRotCOS<TopoCentricRotCOS<BBody, L>> = true;

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopRot   = PosKV   <TopoCentricRotCOS<BBody, L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopRot    = VelKV   <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopRot     = VelV    <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopRot     = AccV    <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopRot   = ForceV  <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopRot  = AngVelV <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopRot  = AngAccV <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopRot  = AngMomV <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopRot    = TorqV   <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopRot     = MoIV    <TopoCentricRotCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopRot = MoIRateV<TopoCentricRotCOS<BBody, L>, B>;

  //=========================================================================//
  // Vector Transforms between the "{Body,Topo}CentricRotCOS"es:             //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // BodyCentric -> TopoCentric:                                             //
  //-------------------------------------------------------------------------//
  // Here the CallER needs to specify the "TCOS" explicitly, which must be some
  // "TopoCentricRotCOS":
  //
  template<typename TCOS, Body  B = Body::UNDEFINED, typename  DQ>
  constexpr Vector3D<DQ,  TCOS, B> ToTopoC
    (Vector3D<DQ, BodyCentricRotCOS<TCOS::BaseBody>, B> const& a_body_c)
  requires(IsTopoCentricRotCOS<TCOS>)
  {
    // Position is changed by parallel translation; all other vectors remain un-
    // changed:
    return
      IsAnyLen<DQ>
      ? Vector3D<DQ, TCOS, B>
        { a_body_c.x() - TCOS::Loc->PosKV().x(),
          a_body_c.y() - TCOS::Loc->PosKV().y(),
          a_body_c.z() - TCOS::Loc->PosKV().z() }
      : Vector3D<DQ, TCOS, B>
        { a_body_c.x(), a_body_c.y(), a_body_c.z() };
  }

  //-------------------------------------------------------------------------//
  // TopoCentric -> BodyCentric:                                             //
  //-------------------------------------------------------------------------//
  // Here the "TopoCentricRotCOS" can be inferred from the argument; the CallER
  // does not need to specify any template params:
  //
  template<Body BBody, Location<BBody> const*  L,
           Body B = Body::UNDEFINED, typename  DQ>
     Vector3D<DQ, BodyCentricRotCOS<BBody>,    B> ToBodyC
    (Vector3D<DQ, TopoCentricRotCOS<BBody, L>, B> const& a_topo_c)
  {
    return
      IsAnyLen<DQ>
      ? PosKVRot<BBody, B>
        { L->PosKV().x() + a_topo_c.x(),
          L->PosKV().y() + a_topo_c.y(),
          L->PosKV().z() + a_topo_c.z() }
      :
        Vector3D<DQ, BodyCentricRotCOS<BBody>, B>
        { a_topo_c.x(), a_topo_c.y(), a_topo_c.z() };
  }

  //=========================================================================//
  // "TopoCentricLaunchCOS" Struct:                                          //
  //=========================================================================//
  // Similar to "TopoCentricRotCOS", but rotated wrt the latter around the Z
  // axis by the angle (counter-clock-wise) (Pi/2-A),  where A is the Launch
  // Azimuth. XXX:  Unfortunately, the latter is dynamic (mission-specific),
  // so it cannot be made into a template param:
  // Often, this COS would be the most convenient one for integration of the
  // Equations of Motion during the Launch phase:
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCentricLaunchCOS
  {
    // XXX: This COS is typically Earth-based, but we allow other Bodies as
    // well, so both TT and TDB are supported:
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    // No objs of this struct are to be created:
    TopoCentricLaunchCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopLaunch    = PosKV   <TopoCentricLaunchCOS<BBody, L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopLaunch    = VelKV   <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopLaunch     = VelV    <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopLaunch     = AccV    <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopLaunch   = ForceV  <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopLaunch  = AngVelV <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopLaunch  = AngAccV <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopLaunch  = AngMomV <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopLaunch    = TorqV   <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopLaunch     = MoIV    <TopoCentricLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopLaunch = MoIRateV<TopoCentricLaunchCOS<BBody, L>, B>;

  //=========================================================================//
  // "AtmCOS" and "VelVAtm":                                                 //
  //=========================================================================//
  // "AtmCOS" is a placeholder used to define "VelVAtm" which is the Velocity
  // Vector relative to the Body's Atmosphere (for various aerodynamic comput-
  // ations).
  // In this case, the "B" Body is certainly an LVSC, so "UNDEFINED": no extra
  // template param is provided:
  //
  template<Body BBody, Location<BBody> const* L>
  struct AtmCOS
  {
    constexpr static  Body BaseBody       = BBody;
    constexpr static  bool HasFixedAxes   = false;
    constexpr static  bool HasFixedOrigin = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;
    static_assert(L != nullptr);
    AtmCOS() = delete; // No objects of this type are to be created
  };

  template<Body BBody, Location<BBody> const* L>
  using VelVAtm  =  VelV<AtmCOS<BBody, L>>;

  // "GetVelVAtm":
  // Velocity vector relative to the Atmosphere ("VelVAtm") can be constructed
  // from the TopoCentric Velocity and TopoCentric Wind Velocity,    in either
  // the "Rot" or "Launch" COSes:
  //
  template<Body BBody, Location<BBody> const* L>
  VelVAtm <BBody, L>   GetVelVAtm
  (
    VelVTopRot   <BBody, L> const& a_vel,
    VelVTopRot   <BBody, L> const& a_wind_vel
  )
  { return a_vel - a_wind_vel; }

  template<Body BBody, Location<BBody> const* L>
  VelVAtm <BBody, L>   GetVelVAtm
  (
    VelVTopLaunch<BBody, L> const& a_vel,
    VelVTopLaunch<BBody, L> const& a_wind_vel
  )
  { return a_vel - a_wind_vel; }
}
// End namespave SpaceBallistics
