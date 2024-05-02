// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/SelenoCentricRotatingCOS.h":         //
//              SelenoCentric Co-Ords System Rotating with the Moon          //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "SelenoCentricRotatingCOS" Class:                                       //
  //=========================================================================//
  // Origin: Center of the Earth ellipsoid
  // Axes  : (X=(0 deg, 0 deg), Y=(90 deg, 0 deg), Z=(*, 90 deg) = NorthPole)
  // NB    : This class stands for itself; no objects of it can be created :
  //
  class SelenoCentricRotatingCOS
  {
    SelenoCentricRotatingCOS() = delete;
  };

  //=========================================================================//
  // Position, Velocity and other Vectors in this COS:                       //
  //=========================================================================//
  using PosVSR    = PosV   <SelenoCentricRotatingCOS>;
  using VelVSR    = VelV   <SelenoCentricRotatingCOS>;
  using AccVSR    = AccV   <SelenoCentricRotatingCOS>;
  using ForceVSR  = ForceV <SelenoCentricRotatingCOS>;
  using AngVelVSR = AngVelV<SelenoCentricRotatingCOS>;
  using AngAccVSR = AngAccV<SelenoCentricRotatingCOS>;
  using AngMomVSR = AngMomV<SelenoCentricRotatingCOS>;
  using TorqVSR   = TorqV  <SelenoCentricRotatingCOS>;
}
// End namespace SpaceBallistics
