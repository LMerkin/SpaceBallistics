// vim:ts=2:et
//===========================================================================//
//                    "Src/PhysForces/DE440T-Data.cpp":                      //
//===========================================================================//
#include "SpaceBallistics/PhysForces/DE440T.h"

namespace SpaceBallistics::DE440T::Bits
{  
  double const Data[NR][ND]   
  {
#   include "SpaceBallistics/PhysForces/DE440T-1650-1750.hpp"
#   include "SpaceBallistics/PhysForces/DE440T-1750-1850.hpp"
#   include "SpaceBallistics/PhysForces/DE440T-1850-1950.hpp"
#   include "SpaceBallistics/PhysForces/DE440T-1950-2050.hpp"
#   include "SpaceBallistics/PhysForces/DE440T-2050-2150.hpp"
  };
}
// End namespace SpaceBallistics::DE440T::Bits