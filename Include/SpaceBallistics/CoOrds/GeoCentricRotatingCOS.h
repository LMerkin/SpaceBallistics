// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/GeoCentricRotatingCOS.h":            //
//              GeoCentric Co-Ords System Rotating with the Earth            //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "GeoCentricRotatingCOS" Class:                                          //
  //=========================================================================//
  // Origin: Center of the Earth ellipsoid
  // Axes  : (X=(0 deg, 0 deg), Y=(90 deg, 0 deg), Z=(*, 90 deg) = NorthPole)
  // NB    : This class stands for itself; no objects of it can be created :
  //
  class GeoCentricRotatingCOS
  {
    GeoCentricRotatingCOS() = delete;
  };

  //=========================================================================//
  // Position, Velocity and other Vectors in this COS:                       //
  //=========================================================================//
  using PosVGR    = PosV   <GeoCentricRotatingCOS>;
  using VelVGR    = VelV   <GeoCentricRotatingCOS>;
  using AccVGR    = AccV   <GeoCentricRotatingCOS>;
  using ForceVGR  = ForceV <GeoCentricRotatingCOS>;
  using AngVelVGR = AngVelV<GeoCentricRotatingCOS>;
  using AngAccVGR = AngAccV<GeoCentricRotatingCOS>;
  using AngMomVGR = AngMomV<GeoCentricRotatingCOS>;
  using TorqVGR   = TorqV  <GeoCentricRotatingCOS>;
}
// End namespace SpaceBallistics
