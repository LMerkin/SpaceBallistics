// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/LVSC/Soyuz-2.1b/Consts.h":              //
//    Over-All Consts for the Mathematical Model of the "Soyuz-2.1b" LV      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics::Soyuz21b_Consts
{
  //=========================================================================//
  // Embedded Co-Ordinate System:                                            //
  //=========================================================================//
  // The following local co-ord system is used. Most X-coords are NEGATIVE but
  // the origin  is not  affected by separation of stages, which is convenient:
  //
  // (*) The OX axis is the main axis of the rocket. The OX positive direction
  //     is towards the HEAD of the LV. The origin O is the LOWER base of the
  //     InterStage (junction plane with Stage3). That is, X >= 0 for
  //     for InterStage, the optional Stage4, Payload Adapter/Dispenser, Pay-
  //     load itself and the Fairing, and X <= 0 for Stages 3, 2, 1.
  // (*) The OY and OZ axes are such that the OXY and OXZ planes pass through
  //     the symmetry axes of the corresp opposite strap-on boosters (Blocks
  //     B, V, G, D -- Stage 1), and OXYZ is a right-oriented co-ords system.
  //     Block B: -Y
  //     Block V: -Z
  //     Block G: +Y
  //     Block D: +Z
  // Source: ArianSpace/StarSem Soyuz CSG User's Manual.

  //=========================================================================//
  // Launch TimeLine:                                                        //
  //=========================================================================//
  using FT = FlightTime;

  // t=0 is the notional LiftOff (Contact Separation) time. Thus, Stage1 and
  // Stage2 engines ignition time is < 0 (the detailed info is in the corresp
  // Stage models):
  constexpr inline FT     LiftOffTime      (0.0_sec);

  // MaxQ occurs at approx 72.0 sec.

  // Stage1 (Blocks B, V, G, D) CutOff and Separation time (earlier versions of
  // the RD-107 engine could burn for up to 140 sec):
  constexpr inline FT     Stage1CutOffTime (118.2_sec);  // Or 117.7 ?
  constexpr inline FT     Stage1SepTime    (118.9_sec);  // Or 118.1 ?

  // Fairing jettisoning time. This is probably NOT a constant, it depends  on
  // mission profile. Some srcs say ~157 sec, ~183 sec, ~228 sec; "Luna-25" was
  // 211.95 sec. So we provide a range instead of a fixed value:
  //
  constexpr inline FT     FairingJetTimeMin(157.0_sec);  // XXX ???
  constexpr inline FT     FairingJetTimeMax(228.0_sec);

  // Stage2 (Block A) separation time. Earlier versions of RD-108 could burn for
  // up to 340 sec:
  constexpr inline FT     Stage2CutOffTime (286.0_sec);  // Or 285?

  // Stage3 ignition occurs at approx "Stage2CutOffTime", PRIOR to Stage2 sepa-
  // ration, immediately at full thrust:
  constexpr inline FT     Stage3IgnTime          = Stage2CutOffTime;
  constexpr inline FT     Stage2SepTime    (287.7_sec);  // 287.6, 287.9?

  // Stage3 aft section jettisoning time;
  // some srcs say (Stage2SepTime + 10 sec, ie ~299 sec):
  constexpr inline FT     Stage3AftJetTime (300.4_sec);

  // Stage3 engine (RD-0124) cut-off time. Some srcs say Stage3 burns for
  // 250..300 sec or even 320 sec  (for the older RD-0110, 240..250 sec);
  // StarSem says 270 sec, "Luna 25" launch was ~274 sec.  XXX: but there are
  // also much shorter burns of 230--240 sec (eg "Prichal" launch); we assume
  // they occur with PARTIAL propellant load in Stage3, which is not supported
  // by our model yet:
  constexpr inline Time   Stage3BurnDur    = 274.0_sec;
  constexpr inline FT     Stage3CutOffTime =
                          Stage3IgnTime + Stage3BurnDur;

  // Payload separation time (from Stage3) -- considered to be the Orbital
  // Insertion time (ArianeSpace/StarSem says 561.9 sec):
  constexpr inline FT     Stage3SepTime    (564.0_sec);
  static_assert(Stage3SepTime > Stage3CutOffTime);

  //=========================================================================//
  // Over-All Dimensions:                                                    //
  //=========================================================================//
  // (So that Stages can, to some extent, be modeled imdependently of each
  // other):
  //
  constexpr inline Len    X0        = 0.0_m; // Top of Stage3
  constexpr inline Len    Stage3Len = 6.745_m;
}
// End namespace SpaceBallistics
