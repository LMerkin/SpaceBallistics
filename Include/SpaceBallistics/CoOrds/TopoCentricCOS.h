// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/TopoCentrticCOS.h":                //
//                      Body-TopoCentric Co-Ord System                      //
//===========================================================================//
#pragma once
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
  template<Body BodyName, Location<BodyName const* L>
  using PosVT     = PosV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using VelVT     = VelV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using AccVT     = AccV   <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using ForceVT   = ForceV <TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using AngVelVT  = AngVelV<TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using AngAccVT  = AngAccV<TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using AngMomVT  = AngMomV<TopoCentricCOS<BodyName, L>>;

  template<Body BodyName, Location<BodyName const* L>
  using TorqVT    = TorqV  <TopoCentricCOS<BodyName, L>>;
}
// End namespave SpaceBallistics
