// vim:ts=2:et
//===========================================================================//
//                             "Soyuz_21b.hpp":                              //
//          Mathematical Model of the "Soyuz-2.1b" Space Launcher            //
//===========================================================================//
#pragma  once
#include "Soyuz_21b.h"
#include "MomentsOfInertia.hpp"
#include "Propellants.hpp"

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
    //-----------------------------------------------------------------------//
    // Launch TimeLine:                                                      //
    //-----------------------------------------------------------------------//
    // Source: AriantSpace/StartSem Soyuz CSG User's Manual, and others.
    // t=0 is the notional LiftOff (Contact Separation) time. Thus, Stage1 and
    // Stage2 engines ignition time is < 0:
    constexpr auto Stages12IgnTime        = Time_sec(-20.0);
    constexpr auto Stages12FullThrustTime = Time_sec_0;

    // MaxQ occurs at approx 72.0 sec.
    //
    // Stage1 (Blocks B, V, G, D) separation time.  XXX: approximately, Stage1
    // engines cut-off occur at the same time  (earlier versions of the RD-107
    // engine could burn for up to 140 sec). More precisely, Stage1 thrust red-
    // uction occurs at 112.0 sec and Stage1 veniers cut-off at 117.7 sec:
    //
    constexpr auto Stage1SepTime          = Time_sec(118.1);

    // Fairing jettisoning time:
    constexpr auto FairingJetTime         = Time_sec(208.4);

    // Stage2 (Block A) separation time. Some srcs indicate 278 sec but this is
    // probably for an earlier version of Soyuz. The RD-108A cut-off time is
    // approx 286 sec. Stage3 ignition occurs at approx that time, just PRIOR to
    // Stage2 separation (but arguably not at full thrust at once). We assume
    // that at the moment of separation, Stage3 is at full thrust.  NB: Earlier
    // versions of RD-108 could burn for up to 340 sec:
    //
    constexpr auto Stage2CutOffTime       = Time_sec(286.0);
    constexpr auto Stage3IgnTime          = Stage2CutOffTime;

    constexpr auto Stage2SepTime          = Time_sec(287.6);
    constexpr auto Stage3FullThrustTime   = Stage2SepTime;

    // Stage3 aft section jettisoning time:
    constexpr auto Stage3AftJetTime       = Time_sec(300.4);

    // Stage3 engine (RD-0124) cut-off time. Some srcs say Stage3 burns for
    // 250..300 sec (the older RD-0110 for 240..250 sec), here we get ~271 sec
    // which is most likely correct:
    constexpr auto Stage3CutOffTime       = Time_sec(558.6);

    // Payload separation time (from Stage3) -- considered to be the Orbital
    // Insertion time:
    constexpr auto Stage3SepTime          = Time_sec(561.9);

    //-----------------------------------------------------------------------//
    // Fairing (Shroud):                                                     //
    //-----------------------------------------------------------------------//
    // XXX: The max diameter of the Fairing may depend on the Payload. However,
    // we currently consider it to be constant:
    constexpr auto FairingDMax            = Len_m(4.11);

    //-----------------------------------------------------------------------//
    // InterStage:                                                           //
    //-----------------------------------------------------------------------//
    constexpr auto InterStageD0           = FairingDMax;
    constexpr auto InterStageD1           = Len_m(2.66);
    constexpr auto InterStageH            = Len_m(0.987);

    //=======================================================================//
    // Stage 3 (Block I):                                                    //
    //=======================================================================//
    // DryMass (incl the jettisonable Aft Section and the RD-0124 engine):
    constexpr auto    Stage3DryMass       = Mass_kg(  2355.0);

    // Mass of the RD-0124 engine (some srcs say 480 kg for an earlier version):
    constexpr auto Stage3EngMass          = Mass_kg(  572.0);

    // Masses of Fuel and Oxidiser. ArianSpace/StarSem says 7600 and 17800 kg,
    // resp., but those figures seem to be incorrect (too much for the actual
    // tank volumes):
    constexpr auto Stage3FuelMass         = Mass_kg( 6650.0);
    constexpr auto Stage3OxidMass         = Mass_kg(16554.0);

    // Over-all Diameter:
    constexpr auto Stage3D                = InterStageD1;

    // Over-all Length of the cylindrical shell:
    constexpr auto Stage3Len              = Len_m(6.745);

    // Fuel Tank, upper (neg-facing) spherical segm: base X-coord and Height:
    constexpr auto Stage3FuelTankUpX0     = Len_m(0.718);
    constexpr auto Stage3FuelTankUpH      = Stage3D/2.0;

    // Fuel Tank, central (cylindrical) section:
    constexpr auto Stage3FuelTankCylX0    = Stage3FuelTankUpX0;
    constexpr auto Stage3FuelTankCylH     = Len_m(0.373);

    // FuelTank,  lower (pos-facing) spherical segm: base X-coord and Height:
    constexpr auto Stage3FuelTankLoX0     = Stage3FuelTankUpX0
                                          + Stage3FuelTankCylH;
    constexpr auto Stage3FuelTankLoH      = Len_m(0.552);

    // Oxidiser Tank, upper (neg-facing) spherical segm: base X-coord and
    // Height:
    constexpr auto Stage3OxidTankUpX0     = Len_m(3.162);
    constexpr auto Stage3OxidTankUpH      = Stage3D/2.0;

    // Oxidiser Tank, central (cylindrical) section:
    constexpr auto Stage3OxidTankCylX0    = Stage3OxidTankUpX0;
    constexpr auto Stage3OxidTankCylH     = Len_m(1.05);

    // Oxidiser Tank, lower (pos-facing) spherical segm: base X-coord and
    // Height:
    constexpr auto Stage3OxidTankLoX0     = Stage3OxidTankUpX0
                                          + Stage3OxidTankCylH;
    constexpr auto Stage3OxidTankLoH      = Stage3D/2.0;

    // Check: The bottom of the FuelTank must not touch the top of the Oxidiser
    // Tank:
    static_assert(Stage3FuelTankLoX0 + Stage3FuelTankLoH <
                  Stage3OxidTankUpX0 - Stage3OxidTankUpH);

    // Also, the bottom of the OxidiserTank must be inside the over-all cylind-
    // rical shell:
    static_assert(Stage3OxidTankLoX0 + Stage3OxidTankLoH < Stage3Len);

    // On the other hand, the top of the FuelTank appears over the over-all
    // cylindrical shell (but inside the InterStage):
    static_assert(IsNeg(Stage3FuelTankUpX0 - Stage3FuelTankUpH)  &&
                        Stage3FuelTankUpX0 - Stage3FuelTankUpH > -InterStageH);

    // Also, compute the volumes of the Fuel and the Oxidiser Tanks and make
    // sure the Propellant masses are within the theorectical maxima:
    constexpr auto Stage3FuelTankVol =
      Volume_SpherSegm(Stage3D, Stage3FuelTankUpH)           +
      Volume_TrCone   (Stage3D, Stage3D, Stage3FuelTankCylH) +
      Volume_SpherSegm(Stage3D, Stage3FuelTankLoH);

    constexpr auto Stage3OxidTankVol =
      Volume_SpherSegm(Stage3D, Stage3OxidTankUpH)           +
      Volume_TrCone   (Stage3D, Stage3D, Stage3OxidTankCylH) +
      Volume_SpherSegm(Stage3D, Stage3OxidTankLoH);

    static_assert(Stage3FuelMass < Stage3FuelTankVol * RG1Dens);
    static_assert(Stage3OxidMass < Stage3OxidTankVol * LOXDens);

    // RD-0124 Engine Characteristics (in Vacuum; SeaLevel is meaningful here):
    constexpr auto Stage3IspVac      = Time_sec(359.0);
    constexpr auto Stage3ThrustVac   = Force   (294.3e3);
    constexpr auto Stage3MassRate    = Stage3ThrustVac / (Stage3IspVac * g0);

    // Then the actual burn time at full thrust should be less than the max
    // theoretical value:
    static_assert(Stage3CutOffTime - Stage3FullThrustTime <
                 (Stage3FuelMass   + Stage3OxidMass) / Stage3MassRate);
  }

  //=========================================================================//
  // Non-Default Ctor:                                                       //
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
	// XXX:
	// (*) All computations here are *APPROXIMATE*, using the thin shell model.
	//     The result is normalised to  the total empty mass of the rocket which
	//     is known reasonably accurately. Engines are modeled as point masses.
	// (*) All consts to be computes are first initialised to NAN, then set to
	//     their actual values via "const_cast":
	//
  template<typename Payload>
  Soyuz_21b<Payload>::Soyuz_21b(Payload const& a_payload)
  : m_payload          (a_payload),
    m_stage3AftMass    (NAN),
    m_stage3EmptyMoI   (NAN)
  {
    using namespace Soyuz_21b_Consts;

    //-----------------------------------------------------------------------//
    // InterStage:                                                           //
    //-----------------------------------------------------------------------//
    // It is a conical shell from FairingDMax to Stage3D, with the OX axis:
    //

    //-----------------------------------------------------------------------//
    // Stage 3 ("Block I"):                                                  //
    //-----------------------------------------------------------------------//
    // Over-all cylindrincal shell with the top at x=0.  This includes the cyl-
    // indrical section of the FuelTank and OxidTank, so make sure thet are not
    // counted twice in the Stage3 MoI:
    auto bigCyl3  =
      MoISV_TrCone_XY(Len_m_0, 0.0, Stage3D, Stage3D, Stage3Len);

    // Stage3 Fuel Tank: UpperSpherSegm (neg-facing), CylMiddle,
    //                   LowerSpherSegm (pos-facing):
    auto fuelUp3  =
      MoISV_SpherSegm_XY
        (false, Stage3FuelTankUpX0, 0.0, Stage3D, Stage3FuelTankUpH);

    auto fuelCyl3 =
      MoISV_TrCone_XY
        (Stage3FuelTankCylX0, 0.0, Stage3D, Stage3D, Stage3FuelTankCylH);

    auto fuelLo3  =
      MoISV_SpherSegm_XY
        (true,  Stage3FuelTankLoX0, 0.0, Stage3D, Stage3FuelTankLoH);

    // Stage3 Oxidiser Tank: Similar to the Fuel Tank:
    auto oxidUp3  =
      MoISV_SpherSegm_XY
        (false, Stage3OxidTankUpX0, 0.0, Stage3D, Stage3OxidTankUpH);

    auto oxidCyl3 =
      MoISV_TrCone_XY
        (Stage3OxidTankCylX0, 0.0, Stage3D, Stage3D, Stage3OxidTankCylH);

    auto oxidLo3  =
      MoISV_SpherSegm_XY
        (true,  Stage3OxidTankLoX0, 0.0, Stage3D, Stage3OxidTankLoH);
	}
}
