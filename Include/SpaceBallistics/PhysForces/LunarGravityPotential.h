// vim:ts=2:et
//===========================================================================//
//           "SpaceBallistics/PhysForces/LunarGravityPotential.h":           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/PhysForces/GravityField.hpp"
#include "SpaceBallistics/CoOrds/SelenoCentricRotatingCOS.h"

namespace SpaceBallistics
{
  // Lunar Gravitational Constant, m^3/sec^2:
  constexpr static GM    KMoon       = GM(4902.8001224453001 * 1e9);

  // Model Equatorial Lunar Radius (m):
  constexpr static Len   ReMoon      = To_Len_m(1738.0_km);

  // The Order of the Gravitational Potential Model:
  constexpr static int   GPOrderMoon = 600;

  // The actual Coeffs are external:
  extern SpherHarmonicCoeffs<SelenoCentricRotatingCOS> const
    GPCoeffsMoon[((GPOrderMoon+1)*(GPOrderMoon+2))/2];
}
// End namespace SpaceBallistics
