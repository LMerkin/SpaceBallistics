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
    Earth,
    Moon,
    Sun,
    Mercury,
    Venus,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune
    // Other Solar System Bodies to come...
  };
}
// End namespace SpaceBallistics
