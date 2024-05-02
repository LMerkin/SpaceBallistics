// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/GeoCentricFixedCOS.h":              //
//               GeoCentric Co-Ords System with ICRF-Fixed Axes              //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "GeoCentricFixedCOS" Class:                                             //
  //=========================================================================//
  // Origin: Center of the Earth ellipsoid.
  // Axes  : X = Equinox(2000.0); XY Plane: Equator(2000.0); Z: by the Equator
  // NB    : This class stands for itself; no objects of it cane be created :
  //
  class GeoCentricFixedCOS
  {
    GeoCentricFixedCOS() = delete;
  };

  //=========================================================================//
  // Position, Velocity and other Vectors in this COS:                       //
  //=========================================================================//
  using PosVGF    = PosV   <GeoCentricFixedCOS>;
  using VelVGF    = VelV   <GeoCentricFixedCOS>;
  using AccVGF    = AccV   <GeoCentricFixedCOS>;
  using ForceVGF  = ForceV <GeoCentricFixedCOS>;
  using AngVelVGF = AngVelV<GeoCentricFixedCOS>;
  using AngAccVGF = AngAccV<GeoCentricFixedCOS>;
  using AngMomVGF = AngMomV<GeoCentricFixedCOS>;
  using TorqVGF   = TorqV  <GeoCentricFixedCOS>;
}
// End namespace SpaceBallistics
