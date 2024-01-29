// vim:ts=2:et
//===========================================================================//
//                              "Soyuz21b_Consts.hpp":                       //
//    Over-All Consts for the Mathematical Model of the "Soyuz-2.1b" LV      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics::Soyuz21b_Consts
{
  //=========================================================================//
  // Launch TimeLine:                                                        //
  //=========================================================================//
  // Source: ArianSpace/StarSem Soyuz CSG User's Manual, and others;
  // t=0 is the notional LiftOff (Contact Separation) time. Thus, Stage1 and
  // Stage2 engines ignition time is < 0:
  constexpr inline Time   Stages12IgnTime        = -15.0_sec;
  constexpr inline Time   Stages12FullThrustTime = 0.0_sec;

  // MaxQ occurs at approx 72.0 sec.

  // Stage1 (Blocks B, V, G, D) separation time.  XXX: approximately, Stage1
  // engines cut-off occur at the same time  (earlier versions of the RD-107
  // engine could burn for up to 140 sec). More precisely, Stage1 thrust red-
  // uction occurs at 112.0 sec and Stage1 veniers cut-off at 117.7 sec:
  //
  constexpr inline Time   Stage1SepTime          = 118.1_sec;

  // Fairing jettisoning time (some srcs say ~183 sec):
  constexpr inline Time   FairingJetTime         = 208.4_sec;

  // Stage2 (Block A) separation time. Some srcs indicate 278 sec but this is
  // probably for an earlier version of Soyuz. The RD-108A cut-off time is
  // approx 286 sec. Stage3 ignition occurs at approx that time, just PRIOR to
  // Stage2 separation (but arguably not at full thrust at once). We assume
  // that at the moment of separation, Stage3 is at full thrust.  NB: Earlier
  // versions of RD-108 could burn for up to 340 sec:
  //
  constexpr inline Time   Stage2CutOffTime       = 286.0_sec;
  constexpr inline Time   Stage3IgnTime          = Stage2CutOffTime;
  constexpr inline Time   Stage3ThrBurnDur       = 1.6_sec;

  constexpr inline Time   Stage2SepTime          = Stage2CutOffTime +
                                                   Stage3ThrBurnDur;
  constexpr inline Time   Stage3FullThrustTime   = Stage2SepTime;

  // Stage3 aft section jettisoning time (some srcs say 297 sec):
  constexpr inline Time   Stage3AftJetTime       = 300.4_sec;

  // Stage3 engine (RD-0124) cut-off time. Some srcs say Stage3 burns for
  // 250..300 sec or even 320 sec  (for the older RD-0110, 240..250 sec),
  // here we get ~271.5 sec at Full Thrust, which is probably about correct:
  constexpr inline Time   Stage3FTBurnDur        = 271.5_sec;
  constexpr inline Time   Stage3ThrottlTime      =
                          Stage3FullThrustTime   + Stage3FTBurnDur;
  constexpr inline Time   Stage3CutOffTime       =
                          Stage3ThrottlTime      + Stage3ThrBurnDur;

  // Payload separation time (from Stage3) -- considered to be the Orbital
  // Insertion time:
  constexpr inline Time   Stage3SepTime          = 561.9_sec;
  static_assert(Stage3SepTime > Stage3CutOffTime);
}
