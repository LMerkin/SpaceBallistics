// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/CoOrds/GeoLocations.hpp":                 //
//                      Constants: Geodetic Locations                        //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz-2*" Launch Sites:                                                //
  //=========================================================================//
  // 1 Launch Pad at Vostochny:
  //
  constexpr static Location<Body::Earth> Vostochny_1S =
    Location<Body::Earth>
    (
      { '+', 128, 20, 05.0 },
      { '+',  51, 53, 03.0 },
      278.0_m
    );

  // 1 active Launch Pad at Baykonur   (the other one, Baykonur_1_5, aka
  // "Gagarin's Start", is currently on conservation and is not suitable
  // for "Soyuz-2*" launches):
  //
  constexpr static Location<Body::Earth> Baykonur_31_6 =
    Location<Body::Earth>
    (
      { '+',  63, 33, 50.0 },
      { '+',  45, 59, 46.0 },
      104.0_m     // Approx
    );

  // 2 active launch Pads as Plesetsk (#3 and #4); #2 was for Molniya-M and has
  // been unused since 2012; and #1 has been demolished in 1999:
  //
  constexpr static Location<Body::Earth> Plesetsk_43_3 =
    Location<Body::Earth>
    (
      { '+',  40, 26, 59.0 },
      { '+',  62, 55, 37.0 },
      100.0_m     // Approx
    );
  constexpr static Location<Body::Earth> Plesetsk_43_4 =
    Location<Body::Earth>
    (
      { '+',  40, 27, 24.0 },
      { '+',  62, 55, 42.0 },
      100.0_m     // Approx
    );
}
// End namespace SpaceBallistics
