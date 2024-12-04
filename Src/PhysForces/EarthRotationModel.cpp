// vim:ts=2:et
//===========================================================================//
//                  "Src/PhysForces/EarthRotationModel.cpp":                 //
//===========================================================================//
#include "SpaceBallistics/PhysForces/EarthRotationModel.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // ITRS -> GCRS Conversions for Standard Vectors:                          //
  //=========================================================================//
  template
  PosKV_GCRS  EarthRotationModel::ToGCRS<LenK> (TT, PosKV_ITRS  const&) const;

  template
  VelKV_GCRS  EarthRotationModel::ToGCRS<VelK> (TT, VelKV_ITRS  const&) const;

  template
  AccV_GCRS   EarthRotationModel::ToGCRS<Acc>  (TT, AccV_ITRS   const&) const;

  template
  ForceV_GCRS EarthRotationModel::ToGCRS<Force>(TT, ForceV_ITRS const&) const;

  //=========================================================================//
  // GCRS -> ITRS Conversions for Standard Vectors:                          //
  //=========================================================================//
  template
  PosKV_ITRS  EarthRotationModel::ToITRS<LenK> (TT, PosKV_GCRS  const&) const;

  template
  VelKV_ITRS  EarthRotationModel::ToITRS<VelK> (TT, VelKV_GCRS  const&) const;

  template
  AccV_ITRS   EarthRotationModel::ToITRS<Acc>  (TT, AccV_GCRS   const&) const;

  template
  ForceV_ITRS EarthRotationModel::ToITRS<Force>(TT, ForceV_GCRS const&) const;
}
// End namespace SpaceBallistics
