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
    // Launch TimeLine:                                                      //
    //=======================================================================//
    // Source: AriantSpace Soyuz CSG User's Manual.
    // t=0 is the notional LiftOff (Contact Separation) time. Thus, Stage1 and
    // Stage2 engines ignition time is < 0. Then:
    //
    // MaxQ occurs at approx 72.0 sec.
    //
    // Stage1 (Blocks B, V, G, D) separation time.  XXX: approximately, Stage1
    // engines cut-off occur at the same time (early versions of the RD-107 en-
    // gine could burn for up to 140 sec). More precisely, Stage1 thrust reduc-
    // tion occurs at 112.0 sec and Stage1 veniers cut-off at 117.7 sec:
    //
    constexpr Time_sec Stage1SepTime    = Time_sec(118.1);

    // Fairing jettisoning time:
    constexpr Time_sec FairingJetTime   = Time_sec(208.4);

    // Stage2 (Block A) separation time. Some srcs indicate 278 sec but this is
    // probably for earlier version of Soyuz. The RD-108A cut-off time is approx
    // 286 sec. Stage3 ignition occurs at approx the latter time, just PRIOR to
    // Stage2 separation (but arguably not at full thrust at once). Early vers-
    // ions of RD-108 could burn for up to 340 sec
    //
    constexpr Time_sec Stage2SepTime    = Time_sec(287.6);

    // Stage3 aft section jettisoning time:
    constexpr Time_sec Stage3AftJetTime = Time_sec(300.4);

    // Stage3 engine (RD-0124) cut-off time. Some srcs say Stage3 burns for
    // 250..300 sec (the older RD-0110 for 240..250 sec), here we get ~271 sec
    // which is most likely correct:
    constexpr Time_sec Stage3CutOffTime = Time_sec(558.6);

    // Payload separation time (from Stage3) -- considered to be the Orbital
    // Insertion time:
    constexpr Time_sec Stage3SepTime    = Time_sec(561.9);

    //=======================================================================//
    // Geometry:                                                             //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Fairing:                                                              //
    //-----------------------------------------------------------------------//
    // XXX: The max diameter of the Fairing may depend on the Payload. However,
    // we currently consider it to be constant:
    constexpr Len_m FairingDMax = Len_m(4.11);

    //-----------------------------------------------------------------------//
    // Stage3:                                                               //
    //-----------------------------------------------------------------------//
    // Diameter:
    constexpr Len_m Stage3D     = Len_m(2.66);

    // Length of the cylindrical shell (not including spherical segments):
    constexpr Len_m Stage3Len   = Len_m(6.745);

    //=======================================================================//
    // Masses:                                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Stage 3 (Block I):                                                    //
    //-----------------------------------------------------------------------//
    // DryMass (incl the jettisonable Aft Section and the RD-0124 engine):
    constexpr auto    Stage3DryMass  = Mass_kg(  2355.0);

    // Mass of the RD-0124 engine (some srcs say 480 for an earlier version):
    constexpr Mass_kg Stage3EngMass  = Mass_kg(  572.0);

    // Masses of Fuel and Oxidiser:
    constexpr Mass_kg Stage3FuelMass = Mass_kg( 7600.0);
    constexpr Mass_kg Stage3OxidMass = Mass_kg(17800.0);

    // Densities of the Fuel (Naftil RG-1) and Oxidiser (LOX):
    constexpr Density Stage3FuelDens = Density( 8330.0);
    constexpr Density Stage3OxidDens = Density(11410.0);
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
    // MoISV sum;

    //-----------------------------------------------------------------------//
    // InterStage:                                                           //
    //-----------------------------------------------------------------------//
    // It is a conical shell from FairingDMax to Stage3D, with the OX axis:
    //
    // sum += MoIS_TrCone(FairingDMax, Stage3D, ...);

    //-----------------------------------------------------------------------//
    // Stage 3 ("Block I"):                                                  //
    //-----------------------------------------------------------------------//

    return MoI(0.0);
	}
}
