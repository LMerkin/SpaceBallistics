// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/Bodies.h":                    //
// Nomenclature of Bodies endowed with Gravitational Field Models and COSes  //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "Body" Enum Class:                                                      //
  //=========================================================================//
  enum class Body: int
  {
    Sun     = 0,
    Mercury = 1,
    Venus   = 2,
    Earth   = 3,
    Mars    = 4,
    Jupiter = 5,
    Saturn  = 6,
    Uranus  = 7,
    Neptune = 8,
    Pluto   = 9,
    // We place the Moon at the end to preserve the classical numbering sequence
    // for the Sun and Planets:
    Moon    = 10
  };
}
// End namespace SpaceBallistics
