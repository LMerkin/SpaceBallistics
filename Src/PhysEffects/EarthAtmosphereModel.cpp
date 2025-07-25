// vim:ts=2:et
//===========================================================================//
//               "Src/PhysEffects/EarthAtmosphereModel.cpp":                 //
//===========================================================================//
// NB: This .CPP file is only required for CLang which cannot evaluate certain
// functions in compile time (as "constexpr"s):
//
#ifdef __clang__
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"

namespace SpaceBallistics::EarthAtmosphereModel
{
  LayerInfo const Layers[7]
  {
# include "SpaceBallistics/PhysEffects/EarthAtmosphereLayers.h"
  };
}
// End namespace SpaceBallistics::EarthAtmosphereModel
#endif
