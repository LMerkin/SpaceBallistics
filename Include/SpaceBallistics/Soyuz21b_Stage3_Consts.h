// vim:ts=2:et
//===========================================================================//
//                        "Soyuz21b_Stage3_Consts.h":                        //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Soyuz21b_Consts.h"
#include "SpaceBallistics/Soyuz21b_Head_Consts.h"

//===========================================================================//
// Consts Related tp Soyuz-2.1b Stage3 ("Block-I"):                          //
//===========================================================================//
namespace SpaceBallistics::Soyuz21b_Stage3_Consts
{
  namespace SHC = Soyuz21b_Head_Consts;
  namespace SC  = Soyuz21b_Consts;

  //-------------------------------------------------------------------------//
  // Geometry:                                                               //
  //-------------------------------------------------------------------------//
  // Stage3 Over-All:
  constexpr inline Len  D               = SHC::InterStageD1;
  constexpr inline Len  H               = 6.745_m;

  // Stage3 Fore Section:
  constexpr inline Len  ForeX0          = 0.0_m;  // THE ORIGIN IS HERE!!!
  constexpr inline Len  ForeH           = 0.718_m;

  // Stage3 Fuel Tank:
  constexpr inline Len  FuelTankMidH    = 0.373_m;
  constexpr inline Len  FuelTankLowH    = 0.552_m;

  // Stage3 Equipment Bay:
  constexpr inline Len  EquipBayH       = 2.071_m;

  // Stage3 Oxidiser Tank:
  constexpr inline Len  OxidTankMidH    = 1.05_m;

  // Stage3 Aft Section:
  constexpr inline Len  AftH            =
    H - (ForeH + FuelTankMidH + EquipBayH + OxidTankMidH);
  static_assert(IsPos(AftH));

  // XXX: The approximate dX of the CoM of the Engine from AftTop:
  constexpr inline Len  EngCoMdX        = 1.0_m;

  //-------------------------------------------------------------------------//
  // Masses:                                                                 //
  //-------------------------------------------------------------------------//
  // EmptyMass (incl the jettisonable Aft Section and the RD-0124 engine):
  // XXX: Some srcs give a much larger EmptyMass: 2597..2710 kg, or a lower
  // one: 2355 kg:
  constexpr inline Mass   EmptyMass     = 2490.0_kg;

  // Mass of the jettisonable Aft Section:
  constexpr inline Mass   AftMass       = 441.0_kg;

  // Mass of the RD-0124 engine:
  // Data in various srcs vary significantly: 450 kg, 460 kg, 520 kg, 572 kg;
  // we assume that the latter figure includes the support structure:
  constexpr inline Mass   EngMass       = 572.0_kg;

  // Masses of Fuel and Oxidiser. ArianSpace/StarSem says 7600 and 17800 kg,
  // resp., but those figures seem to be incorrect  (too much for the actual
  // tank volumes!):
  constexpr inline Mass   FuelMass      = 6650.0_kg;
  constexpr inline Mass   OxidMass      = 16554.0_kg;
  constexpr inline Mass   GasesMass     = 30.0_kg;  // He, air, ...
  constexpr inline Mass   FullMass      = EmptyMass + FuelMass + OxidMass +
                                          GasesMass;
 
  // Remnants of the Fuel and Oxidiser in Stage3 at the engine cut-off time:
  constexpr inline Mass   FuelRem       = 104.0_kg; // StarSem: 98 kg
  constexpr inline Mass   OxidRem       = 167.0_kg; // StarSem: 188..207 kg

  //-------------------------------------------------------------------------//
  // RD-0124 (14D23) Engine Performance:                                     //
  //-------------------------------------------------------------------------//
  // (In Vacuum; SeaLevel is meaningless here). StarSem: Thrust = 297.9e3 N
  constexpr inline Time   IspVac        = 359.0_sec;
  constexpr inline Force  ThrustVac     = Force(294.3e3);
  constexpr inline double Throttling    = 0.5;

  // The actual MassRate (from the RD-0124 data) at Full Thrust.
  // XXX: These figures may not be highly accurate, as it is unclear whether
  // they refer to Naftil or Kerosene. The OxidRate is quoted somewhere as
  // 56.7 kg/sec, which is clearly too low.
  // In different sources, Oxidiser/Fuel Ratio is 2.5..2.6, here 2.50:
  // These figure must be consistent with the BurnDur and the StaticThrust
  // below:
  constexpr inline auto   FuelRateFT    = 23.8_kg / 1.0_sec;
  constexpr inline auto   OxidRateFT    = 59.6_kg / 1.0_sec;
  constexpr inline auto   MassRateFT    = FuelRateFT + OxidRateFT;

  // The MassRate is connected to Specific Impulse and Thrust, but we must
  // take into account that Thrust is a sum of Reactive Force (proportional
  // to MassRate) and the residual inline pressure force at the nozzle exh-
  // aust:
  constexpr inline Force  StaticThrust  = ThrustVac - MassRateFT * IspVac * g0;
  static_assert(IsPos(StaticThrust));
}
