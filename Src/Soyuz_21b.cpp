// vim:ts=2:et
//===========================================================================//
//                              "Soyuz_21b.h":                               //
//          Mathematical Model of the "Soyuz-2.1b" Space Launcher            //
//===========================================================================//
#include "Soyuz_21b.h"
#include "MomentsOfInteria.hpp"

namespace SpaceBallistics
{
	//=========================================================================//
	// "MkEmptyMoI": Moment of Inertia Computation for an Empty Rocket:        //
	//=========================================================================//
	// XXX: All computations are very approximate, using the thin shell model.
	// The results are normalised to the total empty mass of the rocket which
	// is known reasonably accurately. Engines are modeled as point masses:
	//
	MoI_T MkEmptyMoI()
	{
    //-----------------------------------------------------------------------//
    // Stage 3:                                                              //
    //-----------------------------------------------------------------------//
	}
}
