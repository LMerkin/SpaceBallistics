// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/TopoCentrticCOS.h":                //
//                      Earth TopoCentric Co-Ord System                      //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/GeoLocations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TopoCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: A point on the Earth surface given by "L"
  // Axes  : (X=East, Y=North, Z=Zenith)
  // NB    : This class just stands for itself;
  //         no objects of it are to be created:
  //
  template<Location_WGS84 const* L>
  class TopoCentricCOS
  {
    TopoCentricCOS() = delete;
  };

  //=========================================================================//
  // Position, Velocity and other Vectors in this COS:                       //
  //=========================================================================//
  template<Location_WGS84 const* L>
  using PosVET     = PosV   <TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using VelVET     = VelV   <TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using AccVET     = AccV   <TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using ForceVET  `= ForceV <TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using AngVelVET  = AngVelV<TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using AngAccVET  = AngAccV<TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using AngMomVET   = AngMomV<TopoCentricCOS<L>>;

  template<Location_WGS84 const* L>
  using TorqVET    = TorqV  <TopoCentricCOS<L>>;
}
// End namespave SpaceBallistics
