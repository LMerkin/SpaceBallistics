// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/SelenoCentricFixedCOS.h":            //
//              SelenoCentric Co-Ords System with ICRF-Fixed Axes            //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "SelenoCentricFixedCOS" Class:                                          //
  //=========================================================================//
  // Origin: Center of the Earth ellipsoid.
  // Axes  : X = Equinox(2000.0); XY Plane: Equator(2000.0); Z: by the Equator
  // NB    : This class stands for itself; no objects of it cane be created :
  //
  class SelenoCentricFixedCOS
  {
    SelenoCentricFixedCOS() = delete;
  };

  //=========================================================================//
  // Position, Velocity and other Vectors in this COS:                       //
  //=========================================================================//
  using PosVSF    = PosV   <SelenoCentricFixedCOS>;
  using VelVSF    = VelV   <SelenoCentricFixedCOS>;
  using AccVSF    = AccV   <SelenoCentricFixedCOS>;
  using ForceVSF  = ForceV <SelenoCentricFixedCOS>;
  using AngVelVSF = AngVelV<SelenoCentricFixedCOS>;
  using AngAccVSF = AngAccV<SelenoCentricFixedCOS>;
  using AngMomVSF = AngMomV<SelenoCentricFixedCOS>;
  using TorqVSF   = TorqV  <SelenoCentricFixedCOS>;
}
// End namespace SpaceBallistics
