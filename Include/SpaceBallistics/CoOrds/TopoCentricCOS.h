// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/TopoCentrticCOS.h":                //
//                      Body-TopoCentric Co-Ord System                      //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TopoCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: A point on the Body Surface given by "L"
  // Axes  : (X=East, Y=North, Z=Zenith)
  // NB    : This class just stands for itself;
  //         no objects of it are to be created:
  //
  template<Body BodyName, Location<BodyName> const* L>
  class TopoCentricCOS
  {
    static_assert(L != nullptr);
    TopoCentricCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  template<Body BodyName, Location<BodyName> const* L>
  using PosVTop     = PosV    <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName> const* L>
  using VelVTop     = VelV    <TopoCentricCOS<BodyName, L>>;

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
  // "AirCOS" and "VelVAir" Classes:                                         //
  //=========================================================================//
  // "AirCOS" is a placeholder used to define "VelVAir" which is the Air Velo-
  // city Vector (for various aerodynamic computations). XXX: Currently, they
  // are for Earth only, and are unrelated to any particular Location:
  //
  class AirCOS {};
  using VelVAir = VelV<AirCOS>;

  // The Air Velocity can be derived from the TopoCentric Velocity and the
  // (also TopoCentric) Wind Velocity at any Location:
  //
  
}
// End namespave SpaceBallistics
