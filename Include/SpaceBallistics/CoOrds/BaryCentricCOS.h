// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/BaryCentricCOS.h":                  //
//                    The Inertial Co-Ordinate System                        //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "BaryCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: Solar System BaryCenter
  // Axes  : ICRS/ICRF
  // Obviously, using TDB as the assiciated TimeScale:
  //
  class  TDB;

  struct BaryCentricCOS
  {
    using TimeScale  = TDB;

    BaryCentricCOS() = delete;   // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  using PosKVBary = PosKV<BaryCentricCOS>;
  using VelKVBary = VelKV<BaryCentricCOS>;

  // XXX: Currently no need to consider over Vectors in this COS yet...
}
// End namespace SpaceBallistics
