// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/TopoCentrticCOS.h":                //
//                      Body-TopoCentric Co-Ord System                      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/Locations.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "TopoCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin   : A point on the Body Surface given by "L";
  // Axes     : (X=East, Y=North, Z=Zenith);
  // TimeScale: TT for Earth, TDB otherwise;
  // NB       : This class just stands for itself;
  //            no objects of it are to be created:
  //
  class TT;
  class TDB;

  template<Body BodyName, Location<BodyName> const* L>
  struct TopoCentricCOS
  {
    using  TimeScale = std::continional_t<BodyName == Body::Earth ? TT : TDB>;

    static_assert(L != nullptr);
    TopoCentricCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BodyName, Location<BodyName> const* L>
  using PosKVTop    = PosKV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using VelKVTop    = VelKV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AccVTop     = AccV    <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using ForceVTop   = ForceV  <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngVelVTop  = AngVelV <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngAccVTop  = AngAccV <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using AngMomVTop  = AngMomV <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using TorqVTop    = TorqV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIVTop     = MoIV    <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using MoIRateVTop = MoIRateV<TopoCentricCOS<BodyName, L>>;

  //=========================================================================//
  // "AtmCOS" and "VelVAtm" Classes:                                         //
  //=========================================================================//
  // "AtmCOS" is a placeholder used to define "VelVAtm" which is the Velocity
  // Vector relative to the Body's Atmosphere (for various aerodynamic comput-
  // ations):
  //
  template<Body BodyName, Location<BodyName> const* L>
  class AtmCOS
  {
    using  TimeScale = std::continional_t<BodyName == Body::Earth ? TT : TDB>;
    static_assert(L != nullptr);
    AtmCOS() = delete;
  }

  template<Body BodyName, Location<BodyName> const* L>
  using VelVAtm = VelV<AtmCOS<BodyName, L>>;

  // "GetVelVAtm":
  // Velocity vector relative to the Atmosphere ("VelVAtm") can be constructed
  // from the TopoCentric Velocity and TopoCentric Wind Velocity:
  //
  template<Body BodyName, Location<BodyName> const* L>
  VelVAtm <BodyName, L> GetVelVAtm
  (
    VelVTop<BodyName, L> const& a_vel,
    VelVTop<BodyName, L> const& a_wind_vel = VelVTop<BodyName, L>()
  )
  { return a_vel - a_wind_vel; }
}
// End namespave SpaceBallistics
