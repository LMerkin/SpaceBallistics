// vim:ts=2
//===========================================================================//
//                 "SpaceBallistics/LVSC/Propellants.h":                     //
//             Characteristics of Various Rocket Propellants                 //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics::Propellants
{
  // Densities of Naftil (RG-1), Kerosene (T-1) and LOX:
  constexpr auto RG1Dens  = Density( 833.0);	// Naftil
  constexpr auto T1Dens   = Density( 820.0);  // Kerosene; another src: >= 800
  constexpr auto LOxDens  = Density(1141.0);
	constexpr auto H2O2Dens = Density(1420.0);  // 95%
}
// End namespace SpaceBallistics
