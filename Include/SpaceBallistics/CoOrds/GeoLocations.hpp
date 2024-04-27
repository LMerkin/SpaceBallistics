// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/GeoLocations.hpp":                 //
//                      Constants: Geodetic Locations                        //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz-2*" Launch Sites:                                                //
  //=========================================================================//
  // 1 Launch Pad at Vostochny:
  //
  constexpr static Location_WGS84 Vostochny_1S =
    Location_WGS84
    (
      { '+', 128, 20, 05.0 },
      { '+',  51, 53, 03.0 },
      278.0_m
    );

  // 1 active Launch Pad at Baykonur   (the other one, Baykonur_1_5, aka
  // "Gagarin's Start", is currently on conservation and is not suitable
  // for "Soyuz-2*" launches):
  //
  constexpr static Location_WGS84 Baykonur_31_6 =
    Location_WGS84
    (
      { '+',  63, 33, 50.0 },
      { '+',  45, 59, 46.0 },
      104.0_m     // Approx
    );

  // 2 active launch Pads as Plesetsk (#3 and #4); #2 was for Molniya-M and has
  // been unused since 2012 and #1 has been demolished in 1999:
  //
  constexpr static Location_WGS84 Plesetsk_43_3 =
    Location_WGS84
    (
      { '+',  40, 26, 59.0 },
      { '+',  62, 55, 37.0 },
      100.0_m     // Approx
    );
  constexpr static Location_WGS84 Plesetsk_43_4 =
    Location_WGS84
    (
      { '+',  40, 27, 24.0 },
      { '+',  62, 55, 42.0 },
      100.0_m     // Approx
    );
}
