// vim:ts=2:et
//===========================================================================//
//                         "Soyuz2_LaunchPads.hpp":                          //
//                         Locations with Azimuths                           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz2_LaunchPad" Class:                                               //
  //=========================================================================//
  // "Azimuth" is the azimuth (0..360 degs, from the North, Clock-Wise) of the
  // Main Axis of the Pad (ie pointing towards the Far Edge of the Pad):
  //
  class Soyuz2_LaunchPad: public Location
  {
  private:
    Angle   m_mainAxisAz;

  public:
    //-----------------------------------------------------------------------//
    // Ctors:                                                                //
    //-----------------------------------------------------------------------//
  };
}
