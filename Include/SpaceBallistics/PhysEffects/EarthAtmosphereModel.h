// vim:ts=2:et
//===========================================================================//
//           "SpaceBallistics/PhysEffects/EarthAtmosphereModel.h":           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics::EarthAtmosphereModel
{
  //-------------------------------------------------------------------------//
  // Consts:                                                                 //
  //-------------------------------------------------------------------------//
  // Standard Atmospheric Pressure at Sea Level:
  constexpr inline Pressure P0       = Pressure(101325.0);

  // Standard Atmospheric Density  at Sea Level:
  constexpr inline Density  Rho0     = Density(1.225);

  // Specific Gas Constant for Dry Air (J/K):
  constexpr inline auto     RAir     = 287.0528 * Energy(1.0) / 1.0_K;

  // cP/cV Ratio for the Dry Air:
  constexpr inline double   GammaAir = 1.4;
}
// End namespace SpaceBallistics::EarthAtmosphereModel
