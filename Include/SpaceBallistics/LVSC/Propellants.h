// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/LVSC/Propellants.h":                     //
//             Characteristics of Various Rocket Propellants                 //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics::Propellants
{
  // Densities of Naftil (RG-1), Kerosene (T-1) and LOX:
  constexpr inline Density RG1Dens  = Density( 833.0); // Naftil
  constexpr inline Density T1Dens   = Density( 820.0); // Kerosene; alt: >= 800
  constexpr inline Density LOxDens  = Density(1141.0);
  constexpr inline Density H2O2Dens = Density(1345.0); // At 82.5% concentration
  constexpr inline Density LN2Dens  = Density( 806.11);
}
// End namespace SpaceBallistics
