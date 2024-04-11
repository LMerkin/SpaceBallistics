// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/GeoCentricFixedCOS.h":              //
//              GeoCentric Co-Ords System Fixed with the Earth            //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "GeoCentricFixedCOS" Class:                                             //
  //=========================================================================//
  // Origin: Center of the Earth ellipsoid
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
  // We can declare them here because the "GeoCentricROtatingCOS" is a monomor-
  // phic type:
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
