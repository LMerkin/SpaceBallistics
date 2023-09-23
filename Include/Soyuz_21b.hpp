// vim:ts=2:et
//===========================================================================//
//                             "Soyuz_21b.hpp":                              //
//          Mathematical Model of the "Soyuz-2.1b" Space Launcher            //
//===========================================================================//
#pragma  once
#include "Soyuz_21b.h"
#include "MomentsOfInertia.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // Constants:                                                              //
  //=========================================================================//
  // These constants describe the geometry and masses of the Soyuz-2.1b LV. Of
  // course, they could be declared inside the "Soyuz_21b" class. But arguably,
  // they represent implementation details rather than the class spec,  so  we
  // better declare them in a special namespace outside the class:
  //
  namespace Soyuz_21b_Consts
  {
    //=======================================================================//
    // Geometry:                                                             //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Fairing:                                                              //
    //-----------------------------------------------------------------------//
    // XXX: The max diamenter of the Fairing may depend on the Payload. However,
    // we currently consider it to be constant:
    constexpr static Len_m Fairing_DMax = Len_m(4.11);

    //-----------------------------------------------------------------------//
    // Stage3:                                                               //
    //-----------------------------------------------------------------------//
    // Diameter:
    constexpr static Len_m Stage3_D     = Len_m(2.66);

    // Length of the cylindrical shell (not including spherical segments):
    constexpr static Len_m Stage3_Len   = Len_m(6.745);

    //=======================================================================//
    // Masses:                                                               //
    //=======================================================================//
  }

  //=========================================================================//
  // Non-Default Ctor:                                                       //
  //=========================================================================//
  template<typename Payload>
  Soyuz_21b<Payload>::Soyuz_21b(Payload const& a_payload)
  : m_payload (a_payload),
    m_emptyMoI(EmptyMoI())
  {}

	//=========================================================================//
	// "EmptyMoI": Moment of Inertia Computation for the Empty Rocket:         //
	//=========================================================================//
  // For MoI computations, the following local co-ord system is used for conve-
  // nience (so that most X-coords are positive):
  // (*) The OX axis is the main axis of the rocket. The positive direction is
  //     towards the TAIL (REAL of the LV). The origin O  is at the connection
  //     plane between the InterStage and Stage3. That is, X < 0 for the rocket
  //     head (InterStage, optional Stage4, Payload and Fairing), and X > 0 for
  //     the most part of Stage3, for Stage2 and Stage1.
  // (*) The OY and OZ axes are such that the OXY and OXZ planes pass through
  //     the symmetry axes of the corresp opposite strap-on boosters (Stage 1),
  //     and OXYZ is a right-oriented co-ords system.
  //
	// XXX: All computations here are *APPROXIMATE*, using the thin shell model.
	// The result is normalised to  the total empty mass of the rocket which is
	// known reasonably accurately. Engines are modeled as point masses.
	//
	template<typename Payload>
	MoI Soyuz_21b<Payload>::EmptyMoI() const
	{
    using namespace Soyuz_21b_Consts;
    //-----------------------------------------------------------------------//
    // Shells:                                                               //
    //-----------------------------------------------------------------------//
    // Here we accumulate MoI per Surface Density and the Surface Area of all
    // Shells (initially 0s):
    MoIS sum;

    //-----------------------------------------------------------------------//
    // InterStage:                                                           //
    //-----------------------------------------------------------------------//
    // It is a conical shell from FairingDMax to Stage3D, with the OX axis:
    //
    // sum += MoIS_TrCone(FairingDMax, Stage3D

    //-----------------------------------------------------------------------//
    // Stage 3 ("Block I"):                                                  //
    //-----------------------------------------------------------------------//

    return MoI(0.0);
	}
}
