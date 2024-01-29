// vim:ts=2
//===========================================================================//
//                             "Propellants.h":                              //
//             Characteristics of Various Rocket Propellants                 //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  // Densities of Naftil (RG-1), Kerosene (T-1) and LOX:
  constexpr auto RG1Dens = Density( 833.0);
  constexpr auto T1Dens  = Density( 820.0);  // Another src: >= 800
  constexpr auto LOxDens = Density(1141.0);
}
