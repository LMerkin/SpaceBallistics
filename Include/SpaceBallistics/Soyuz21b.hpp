// vim:ts=2:et
//===========================================================================//
//                              "Soyuz21b.hpp":                              //
//          Mathematical Model of the "Soyuz-2.1b" Launch Vehicle            //
//===========================================================================//
// Geometrical, Mass and Engine Constant Params.
//
// For MoI computations, the following local co-ord system is used for con-
// venience (so that most X-coords are positive):
// (*) The OX axis is the main axis of the rocket. The OX positive direction
//     is towards the TAIL (REAL of the LV). The origin O is the upper base
//     of the InterStage (junction plane with the Fairing). That is, X < 0
//     for the optional Stage4, payload adapter/dispenser,  payload itself
//     and the fairing, and X >= 0 for the Integstage and Stages 3, 2, 1.
// (*) The OY and OZ axes are such that the OXY and OXZ planes pass through
//     the symmetry axes of the corresp opposite strap-on boosters (Blocks
//     B, V, G, D -- Stage 1), and OXYZ is a right-oriented co-ords system.
// XXX:
// (*) All computations based on thses params are *APPROXIMATE*, using the
//     thin shell model normalised to the empty masses of the rocket stages
//     which are known reasonably accurately...
//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/Propellants.h"
#include "SpaceBallistics/ConstrElement.hpp"

namespace SpaceBallistics::Soyuz21b_Consts
{
  using CE = ConstrElement;

  //=========================================================================//
  // Launch TimeLine:                                                        //
  //=========================================================================//
  // Source: ArianSpace/StarSem Soyuz CSG User's Manual, and others;
  // t=0 is the notional LiftOff (Contact Separation) time. Thus, Stage1 and
  // Stage2 engines ignition time is < 0:
  constexpr Time   Stages12IgnTime        = -15.0_sec;
  constexpr Time   Stages12FullThrustTime = 0.0_sec;

  // MaxQ occurs at approx 72.0 sec.

  // Stage1 (Blocks B, V, G, D) separation time.  XXX: approximately, Stage1
  // engines cut-off occur at the same time  (earlier versions of the RD-107
  // engine could burn for up to 140 sec). More precisely, Stage1 thrust red-
  // uction occurs at 112.0 sec and Stage1 veniers cut-off at 117.7 sec:
  //
  constexpr Time   Stage1SepTime          = 118.1_sec;

  // Fairing jettisoning time (some srcs say ~183 sec):
  constexpr Time   FairingJetTime         = 208.4_sec;

  // Stage2 (Block A) separation time. Some srcs indicate 278 sec but this is
  // probably for an earlier version of Soyuz. The RD-108A cut-off time is
  // approx 286 sec. Stage3 ignition occurs at approx that time, just PRIOR to
  // Stage2 separation (but arguably not at full thrust at once). We assume
  // that at the moment of separation, Stage3 is at full thrust.  NB: Earlier
  // versions of RD-108 could burn for up to 340 sec:
  //
  constexpr Time   Stage2CutOffTime       = 286.0_sec;
  constexpr Time   Stage3IgnTime          = Stage2CutOffTime;

  constexpr Time   Stage2SepTime          = 287.6_sec;
  constexpr Time   Stage3FullThrustTime   = Stage2SepTime;

  // Stage3 aft section jettisoning time (some srcs say 297 sec):
  constexpr Time   Stage3AftJetTime       = 300.4_sec;

  // Stage3 engine (RD-0124) cut-off time. Some srcs say Stage3 burns for
  // 250..300 sec or even 320 sec  (for the older RD-0110, 240..250 sec),
  // here we get ~271 sec which is probably correct:
  constexpr Time   Stage3CutOffTime       = 558.6_sec;

  // Payload separation time (from Stage3) -- considered to be the Orbital
  // Insertion time:
  constexpr Time   Stage3SepTime          = 561.9_sec;

  //=========================================================================//
  // Fairing (Shroud):                                                       //
  //=========================================================================//
  // XXX: The max diameter of the Fairing may depend on the Payload. However,
  // we currently consider it to be constant. We assume the LARGEST Fairing:
  //
  constexpr Len    FairingDMax      = 4.11_m;   // Small fairing: 2.7 m

  // Radius of the Fairing' sperical top:
  constexpr Len    FairingSpherR    = 1.18_m;

  // Half-angle of the Fairing's conical part: 13.5 deg:
  constexpr double FairingHalfAngle = 13.5 * M_PI / 180.0;

  // Over-all Fairing length, incl the spherical (top) conical (middle) and
  // cylindrical (lower) part:
  constexpr Len    FairingTotalLen  = 11.433_m; // Small fairing: 8.34 m

  // Over-all Fairing mass:
  constexpr Mass   FairingMass      = 1700.0_kg;

  //=========================================================================//
  // InterStage:                                                             //
  //=========================================================================//
  // InterStage is a truncated conical shell. The upper diameter is that of
  // the Fairing, the lower one -- of Stage3:
  constexpr Len    InterStageX0     = 0.0_m;    // The origin of X
  constexpr Len    InterStageD0     = FairingDMax;
  constexpr Len    InterStageD1     = 2.66_m;
  constexpr Len    InterStageH      = 1.20_m;   // Small fairing: 0.987 m

  constexpr Mass   InterStageMass   = 400.0_kg;

  //=========================================================================//
  // Stage 3 (Block I):                                                      //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Masses:                                                                 //
  //-------------------------------------------------------------------------//
  // EmptyMass (incl the jettisonable Aft Section and the RD-0124 engine):
  // XXX: Some srcs give a much larger EmptyMass: 2597..2710 kg; maybe the
  // latter includes the InterStage?
  constexpr Mass   Stage3EmptyMass     = 2355.0_kg;

  // Mass of the RD-0124 engine (some srcs say 480 kg for an earlier version):
  constexpr Mass   Stage3EngMass       = 572.0_kg;

  // Masses of Fuel and Oxidiser. ArianSpace/StarSem says 7600 and 17800 kg,
  // resp., but those figures seem to be incorrect  (too much for the actual
  // tank volumes!):
  constexpr Mass   Stage3FuelMass      = 6650.0_kg;
  constexpr Mass   Stage3OxidMass      = 16554.0_kg;
  constexpr Mass   Stage3HeMass        = 27.0_kg;

  // Remnants of the Fuel and Oxidiser in Stage3 at the engine cut-off time:
  constexpr Mass   Stage3FuelRem       = 104.0_kg; // A: 98 kg
  constexpr Mass   Stage3OxidRem       = 167.0_kg; // A: 188..207 kg

  //-------------------------------------------------------------------------//
  // Geometry:                                                               //
  //-------------------------------------------------------------------------//
  // Over-all Diameter:
  constexpr Len    Stage3D             = InterStageD1;

  // Over-all Base and Length of the cylindrical shell:
  constexpr Len    Stage3X0            = InterStageX0 + InterStageH;
  constexpr Len    Stage3H             = 6.745_m;

  // Fuel Tank, upper (neg-facing) spherical segm: base X-coord and Height:
  constexpr Len    Stage3FuelTankUpX0  = Stage3X0 + 0.718_m;
  constexpr Len    Stage3FuelTankUpH   = Stage3D/2.0;

  // Fuel Tank, central (cylindrical) section:
  constexpr Len    Stage3FuelTankCylX0 = Stage3FuelTankUpX0;
  constexpr Len    Stage3FuelTankCylH  = 0.373_m;

  // FuelTank,  lower (pos-facing) spherical segm: base X-coord and Height:
  constexpr Len    Stage3FuelTankLoX0  = Stage3FuelTankCylX0
                                       + Stage3FuelTankCylH;
  constexpr Len    Stage3FuelTankLoH   = 0.552_m;

  // Oxidiser Tank, upper (neg-facing) spherical segm: base X-coord and
  // Height:
  constexpr Len    Stage3OxidTankUpX0  = Stage3X0 + 3.162_m;
  constexpr Len    Stage3OxidTankUpH   = Stage3D/2.0;

  // Oxidiser Tank, central (cylindrical) section:
  constexpr Len    Stage3OxidTankCylX0 = Stage3OxidTankUpX0;
  constexpr Len    Stage3OxidTankCylH  = 1.05_m;

  // Oxidiser Tank, lower (pos-facing) spherical segm: base X-coord and
  // Height:
  constexpr Len    Stage3OxidTankLoX0  = Stage3OxidTankCylX0
                                       + Stage3OxidTankCylH;
  constexpr Len    Stage3OxidTankLoH   = Stage3D/2.0;

  // The Jettisonable Aft Section:
  constexpr Len    Stage3AftX0         = Stage3OxidTankLoX0;
  constexpr Len    Stage3AftH          = Stage3X0 + Stage3H - Stage3AftX0;

  // FIXME: The approximate X-coord of the CoM of the Engine:
  constexpr Len    Stage3EngX0         = Stage3AftX0 + 1.0_m;

  //-------------------------------------------------------------------------//
  // RD-0124 (14D23) Engine Performance:                                     //
  //-------------------------------------------------------------------------//
  // (In Vacuum; SeaLevel is meaningless here):
  constexpr Time   Stage3IspVac    = 359.0_sec;
  constexpr Force  Stage3ThrustVac = Force(294.3e3);
  // (StarSems says Thrust = 297.9e3 N)

  // The actual MassRate:
  constexpr auto   Stage3MassRate  = 80.6_kg / 1.0_sec;

  // The MassRate is connected to Specific Impulse and Thrust, but we must take
  // into account that Thrust is a sum of Reactive Force (proportional to Mass-
  // Rate) and the residual static pressure force at the nozzle exhaust:
  constexpr auto Stage3StaticThrust =
    Stage3ThrustVac - Stage3MassRate * Stage3IspVac * g0;
  static_assert(IsPos(Stage3StaticThrust));

  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  // The bottom of the FuelTank must not touch the top of the Oxidiser Tank:
  static_assert(Stage3FuelTankLoX0 + Stage3FuelTankLoH <
                Stage3OxidTankUpX0 - Stage3OxidTankUpH);

  // On the other hand, the top of the FuelTank appears over the over-all
  // cylindrical shell (but inside the InterStage):
  static_assert(Stage3FuelTankUpX0 - Stage3FuelTankUpH < Stage3X0   &&
                Stage3FuelTankUpX0 - Stage3FuelTankUpH > InterStageX0);

  // Also, the bottom of the OxidiserTank must be well inside the over-all
  // cylindrical shell:
  static_assert(Stage3OxidTankLoX0 + Stage3OxidTankLoH <
                Stage3X0           + Stage3H);
  // Or equivalently:
  static_assert(Stage3OxidTankLoH  < Stage3AftH);

  // Also, compute the volumes of the Fuel and the Oxidiser Tanks and make
  // sure the Propellant masses are within the theorectical maxima:
  //
  constexpr auto Stage3FuelTankVol =
    CE::Volume_SpherSegm(Stage3D, Stage3FuelTankUpH)           +
    CE::Volume_TrCone   (Stage3D, Stage3D, Stage3FuelTankCylH) +
    CE::Volume_SpherSegm(Stage3D, Stage3FuelTankLoH);

  constexpr auto Stage3OxidTankVol =
    CE::Volume_SpherSegm(Stage3D, Stage3OxidTankUpH)           +
    CE::Volume_TrCone   (Stage3D, Stage3D, Stage3OxidTankCylH) +
    CE::Volume_SpherSegm(Stage3D, Stage3OxidTankLoH);

  static_assert(Stage3FuelMass < Stage3FuelTankVol * RG1Dens);
  static_assert(Stage3OxidMass < Stage3OxidTankVol * LOXDens);

  constexpr auto Stage3SpentMass =
    (Stage3FuelMass + Stage3OxidMass) -
    (Stage3FuelRem  + Stage3OxidRem);

  static_assert(Stage3CutOffTime - Stage3FullThrustTime <
                Stage3SpentMass  / Stage3MassRate);
}
