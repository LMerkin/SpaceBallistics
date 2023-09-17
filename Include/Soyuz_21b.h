// vim:ts=2:et
//===========================================================================//
//                              "Soyuz_21b.h":                               //
//          Mathematical Model of the "Soyuz-2.1b" Space Launcher            //
//===========================================================================//
#pragma  once
#include "Types.h"

namespace SpaceBallistics
{
	//=========================================================================//
	// "Soyuz_21b" Class:																											 //
	//=========================================================================//
	// "Soyuz_21b" considered to be a 3-Stage rocket. The class is parametrised
	// by the Payload (which may in particular include an Upper (4th) Stage, eg
	// "Fregat").
	// The embedded co-ord system Oxyz:
	// (*) The Ox axis is the main axis of the rocket. The positive direction is
	//     towards the nose. The origin O is in the center of the Payload Adapter
	//     Ring.
	// (*) The Oy and Oz axes are such that the Oxy and Oxz planes pass through
	//     the symmetry axes of the corresp opposite strap-on boosters (Stage 1),
	//     and Oxyz is a right-oriented co-ords system.
	//
	template<typename Payload>
  class Soyuz_21b
  {
  private:
    //-----------------------------------------------------------------------//
    // Moment of Inertia Computation for an Empty Rocket:                    //
    //-----------------------------------------------------------------------//
    static MoI_T MkEmptyMoI();
  };
}
