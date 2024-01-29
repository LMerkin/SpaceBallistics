// vim:ts=2:et
//===========================================================================//
//                          "Soyuz21b_Head_Consts.h":                        //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"

//===========================================================================//
// Consts Related tp Soyuz-2.1b InterStage and Fairing, in 2 Variants:       //
//===========================================================================//
namespace SpaceBallistics::Soyuz21b_Head_Consts
{
  //-------------------------------------------------------------------------//
  // Fairing Geometry and Masses: Large and Small:                           //
  //-------------------------------------------------------------------------//
  // Max Diameter:
  constexpr inline Len    FairingLargeMaxD    = 4.11_m;
  constexpr inline Len    FairingSmallMaxD    = 2.70_m;

  // Over-All Length:
  constexpr inline Len    FairingLargeH       = 11.433_m;
  constexpr inline Len    FairingSmallH       =  8.340_m;

  // Detailed geometry of the Large Fairing:
  constexpr inline Len    FairingLargeTopR    = 1.180_m;
  constexpr inline double FairingLargeTopAng  = 13.5 * Pi<double> / 180.0;

  // Detailed geometry of the Small Fairing:
  constexpr inline double FairingSmallTopAng  = 30.0 * Pi<double> / 180.0;

  // NB: ArianeSpace/StarSem says 1700 kg for Large Fairing:
  constexpr inline Mass   FairingLargeMass    = 1550.0_kg;
  constexpr inline Mass   FairingSmallMass    = 1000.0_kg;

  //-------------------------------------------------------------------------//
  // InterStage Geometry and Masses: Large and Small:                        //
  //-------------------------------------------------------------------------//
  constexpr inline Len  InterStageLargeD0     = FairingLargeMaxD;
  constexpr inline Len  InterStageSmallD0     = FairingSmallMaxD;
  constexpr inline Len  InterStageLargeH      = 1.20_m;
  constexpr inline Len  InterStageSmallH      = 0.987_m;
  // The lower diameter is the same for both types of Fairing/InterStage, as it
  // coincides with the Stage3 diameter:
  constexpr inline Len  InterStageD1          = 2.66_m;

  // NB: ArianeSpace/StarSem says 400 kg for Large InterStage:
  constexpr inline Mass InterStageLargeMass   = 450.0_kg;
  constexpr inline Mass InterStageSmallMass   = 220.0_kg;
}
