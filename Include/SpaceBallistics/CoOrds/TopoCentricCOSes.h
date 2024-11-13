// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/TopoCentrticCOSes.h":               //
//                      Body-TopoCentric Co-Ord Systems                      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
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
  // ICRS), but the Origin is a Location ("L") on the Body, not the Body Center.
  // Suitable for computation of TopoCentric Ephemerides in astronomical appls.
  // The Origin of this COS is rotating (and otherwise moving) in the inertial
  // space, but the Axes are fixed:
  //
  template<Body BodyName, Location<BodyName> const* L>
  struct TopoCentricEqFixCOS
  {
    constexpr static  Body BaseBody       = BodyName;
    constexpr static  bool HasFixedAxes   = true;
    constexpr static  bool HasFixedOrigin = false;
    using TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCentricEqFixCOS() = delete; // No objs of this struct are to be created
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors in this COS:                              //
  //-------------------------------------------------------------------------//
  // We probably do not need any other DQs here:
  //
  template<Body BodyName, Location<BodyName> const* L>
  using PosKVTopEqFix = PosKV<TopoCentricEqFixCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using VelKVTopEqFix = VelKV<TopoCentricEqFixCOS<BodyName, L>>;

  //=========================================================================//
  // "TopoCentricRotCOS" Struct:                                             //
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
  template<Body BodyName, Location<BodyName> const* L>
  struct TopoCentricRotCOS
  {
    constexpr static  Body BaseBody       = BodyName;
    constexpr static  bool HasFixedAxes   = false;
    constexpr static  bool HasFixedOrigin = false;
    using TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCentricRotCOS() = delete;  // No objs of this struct are to be created
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BodyName, Location<BodyName> const* L>
  using PosKVTopRot    = PosKV   <TopoCentricRotCOS<BodyName, L>>;

  // "VelK" and "Vel":
  template<Body BodyName, Location<BodyName> const* L>
  using VelKVTopRot    = VelKV   <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using VelVTopRot     = VelV    <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AccVTopRot     = AccV    <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using ForceVTopRot   = ForceV  <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngVelVTopRot  = AngVelV <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngAccVTopRot  = AngAccV <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngMomVTopRot  = AngMomV <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using TorqVTopRot    = TorqV   <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIVTopRot     = MoIV    <TopoCentricRotCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIRateVTopRot = MoIRateV<TopoCentricRotCOS<BodyName, L>>;

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
  template<Body BodyName, Location<BodyName> const* L>
  struct TopoCentricLaunchCOS
  {
    // XXX: This COS is typically Earth-based, but we allow other Bodies as
    // well, so both TT and TDB are supported:
    constexpr static  Body BaseBody       = BodyName;
    constexpr static  bool HasFixedAxes   = false;
    constexpr static  bool HasFixedOrigin = false;
    using TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    // No objs of this struct are to be created:
    TopoCentricLaunchCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BodyName, Location<BodyName> const* L>
  using PosKVTopLaunch    = PosKV   <TopoCentricLaunchCOS<BodyName, L>>;

  // "VelK" and "Vel":
  template<Body BodyName, Location<BodyName> const* L>
  using VelKVTopLaunch    = VelKV   <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using VelVTopLaunch     = VelV    <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AccVTopLaunch     = AccV    <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using ForceVTopLaunch   = ForceV  <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngVelVTopLaunch  = AngVelV <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngAccVTopLaunch  = AngAccV <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngMomVTopLaunch  = AngMomV <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using TorqVTopLaunch    = TorqV   <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIVTopLaunch     = MoIV    <TopoCentricLaunchCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIRateVTopLaunch = MoIRateV<TopoCentricLaunchCOS<BodyName, L>>;

  //=========================================================================//
  // "AtmCOS" and "VelVAtm":                                                 //
  //=========================================================================//
  // "AtmCOS" is a placeholder used to define "VelVAtm" which is the Velocity
  // Vector relative to the Body's Atmosphere (for various aerodynamic comput-
  // ations):
  //
  template<Body BodyName, Location<BodyName> const* L>
  struct AtmCOS
  {
    constexpr static  Body BaseBody       = BodyName;
    constexpr static  bool HasFixedAxes   = false;
    constexpr static  bool HasFixedOrigin = false;
    using TimeScale = std::conditional_t<BodyName == Body::Earth, TT, TDB>;
    static_assert(L != nullptr);
    AtmCOS() = delete; // No objects of this type are to be created
  };

  template<Body BodyName, Location<BodyName> const* L>
  using VelVAtm = VelV<AtmCOS<BodyName, L>>;

  // "GetVelVAtm":
  // Velocity vector relative to the Atmosphere ("VelVAtm") can be constructed
  // from the TopoCentric Velocity and TopoCentric Wind Velocity,    in either
  // the "Rot" or "Launch" COSes:
  //
  template<Body BodyName, Location<BodyName> const* L>
  VelVAtm <BodyName, L> GetVelVAtm
  (
    VelVTopRot   <BodyName, L> const& a_vel,
    VelVTopRot   <BodyName, L> const& a_wind_vel
  )
  { return a_vel - a_wind_vel; }

  template<Body BodyName, Location<BodyName> const* L>
  VelVAtm <BodyName, L> GetVelVAtm
  (
    VelVTopLaunch<BodyName, L> const& a_vel,
    VelVTopLaunch<BodyName, L> const& a_wind_vel
  )
  { return a_vel - a_wind_vel; }
}
// End namespave SpaceBallistics
