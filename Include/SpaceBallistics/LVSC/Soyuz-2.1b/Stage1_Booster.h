// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h":           //
//      Mathematical Model of the "Soyuz-2.1b" StrapOn Booster of Stage1     //
//                        (any of "Blocks B,V,G,D")                          //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/TrConeSpherSegm.hpp"
#include "SpaceBallistics/ME/ToricSegms.hpp"
#include "SpaceBallistics/LVSC/Propellants.h"
#include "SpaceBallistics/LVSC/StageDynParams.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Consts.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC  = Soyuz21b_Consts;

  // All "MechElements" are instantiated with "LVSC::Soyuz21b":
  using     ME  = MechElement<LVSC::Soyuz21b>;
  using     PM  = PointMass  <LVSC::Soyuz21b>;
  using     TrC = TrCone     <LVSC::Soyuz21b>;
  using     SpS = SpherSegm  <LVSC::Soyuz21b>;
  using     Tor = ToricSegm  <LVSC::Soyuz21b>;
  using     S2  = Soyuz21b_Stage2;

  //=========================================================================//
  // "Soyuz21b_Stage1_Booster" Class:                                        //
  //=========================================================================//
  template<char Block>
  class Soyuz21b_Stage1_Booster
  {
  private:
    //=======================================================================//
    // No objects can be created of this class:                              //
    //=======================================================================//
    Soyuz21b_Stage1_Booster() = delete;

  public:
    //=======================================================================//
    // Geometry:                                                             //
    //=======================================================================//
    static_assert
      (Block == 'B' || Block == 'V' || Block == 'G' || Block == 'D');
    // Block B: -Y, Psi=Pi
    // Block V: -Z, Psi=3*Pi/2
    // Block G: +Y, Psi=0
    // Block D: +Z, Psi=Pi/2
    //
    // The Psi angle in the OYZ plane (see "RotationShell" for the definition):
    // Psi = Pi/2 * I:
    constexpr static double   CosPsi =
      (Block == 'B') ? -1.0 : (Block == 'G') ? 1.0 : 0.0;
    constexpr static double   SinPsi =
      (Block == 'V') ? -1.0 : (Block == 'D') ? 1.0 : 0.0;

    //-----------------------------------------------------------------------//
    // The Top:                                                              //
    //-----------------------------------------------------------------------//
    // The X-coord of the Booster top: Relative to MaxD of Stage2.  The "TopX"
    // co-ord is adjusted to align the LoX points of Stage2 and Stage1. "AttH"
    // is the distance below the "MaxD" belt of Stage2 where Stage1 is attached
    // to Stage2, "TopR" is the corresp distance from the OX axis,  and  hence
    // the "TopY" and "TopZ" co-ords:
    //
    static_assert(S2::OxidTankUp.GetLow()[0] == S2::OxidTankLow.GetUp()[0]);
    constexpr static Len      AttH   = 0.08_m;
    constexpr static Len      TopX   = S2::OxidTankUp.GetLow()[0] - AttH;
    constexpr static Len      TopR   =
      S2::MaxD + (S2::MinD - S2::MaxD) * double(AttH / S2::MidTrCH);
    constexpr static Len      TopY   = TopR * CosPsi;
    constexpr static Len      TopZ   = TopR * SinPsi;

    // XXX: The following geometry of is somewhat approximate, it does not take
    // into account the full complicated shape of the Booster Block...
    // The TopCone (assumed to be a full cone, not truncated!):
    constexpr static Len      TopConeH     = 4.015_m;
    constexpr static Len      TopConeD     = 1.350_m;

    // The Top-most part of the TopCone is by itself empty, and is capped with a
    // ball joint used for attaching the StrapOn Booster to the core Stage2:
    constexpr static double   TopTopConeF  = 0.3;
    constexpr static Len      TopTopConeH  = TopTopConeF * TopConeH;
    constexpr static Len      TopTopConeD  = TopTopConeF * TopConeD;

    //-----------------------------------------------------------------------//
    // OxidTank:                                                             //
    //-----------------------------------------------------------------------//
    // At the bottom of the TopTopCone, there is a small SpherSegm which is the
    // Top of the OxidTank:
    constexpr static Len      OxidTankTopH = 0.04 * TopConeH;
    constexpr static Len      OxidTankTopD = TopTopConeD;

    // The low part of the "TopTrCone" contains the Up (TrCone) part of the Oxid
    // Tank:
    constexpr static Len      OxidTankUpH  = TopConeH - TopTopConeH;
    constexpr static Len      OxidTankUpD  = TopConeD;
 
    // Below is the Main ("Mid") part of the OxidTank, which is another TrCone.
    // Its diameter increases from "OxidTankUpD" to "OxidTankMidD":
    constexpr static Len      MaxD         = 2.680_m;
    constexpr static Len      MainTrCH     = 11.565_m; // XXX: Manually Adjusted
    constexpr static Len      OxidTankMidH = 8.25_m;   // ditto
    constexpr static Len      OxidTankMidD =
      OxidTankUpD  + (MaxD - OxidTankUpD)  * double(OxidTankMidH / MainTrCH);

    // Below is the InterTank section which is empty by itself and serves as an
    // enclosure for the OxidTankBtm and FuelTankTop SpherSegms:
    constexpr static Len      InterTankH   = 1.075_m;
    constexpr static Len      InterTankUpD = OxidTankMidD;
    constexpr static Len      InterTankLoD =
      InterTankUpD + (MaxD - OxidTankUpD)  * double(InterTankH   / MainTrCH);

    // OxidTankBtm:
    constexpr static Len      OxidTankBtmH = 0.48_m;
    constexpr static Len      OxidTankBtmD = OxidTankMidD;

    //-----------------------------------------------------------------------//
    // FuelTank:                                                             //
    //-----------------------------------------------------------------------//
    // FuelTankTop:
    constexpr static Len      FuelTankTopH = OxidTankBtmH;
    static_assert(IsPos(InterTankH - OxidTankBtmH - FuelTankTopH));
    constexpr static Len      FuelTankTopD = InterTankLoD;

    // FuelTankMid: Its up diameter is "FuelTankTopD" and the low diameter is
    // "MaxD":
    constexpr static Len      FuelTankMidH =
      MainTrCH - OxidTankMidH - InterTankH;
    static_assert(IsPos(FuelTankMidH));

    // Then the cylindrical enclosure of the FuelTankBtm, H2O2Tank, LiqN2 Tank
    // and the Engine, of the "MaxD" diameter. XXX: We ignore the fact that it
    // is not strictly cylindrical in its lowest (Engine) part. AeroFin is also
    // included into the mass of "TailCylEnc":
    constexpr static Len      TailCylEnclH = 4.020_m;

    // FuelTankBtm:
    constexpr static Len      FuelTankBtmH = OxidTankBtmH;

    //-----------------------------------------------------------------------//
    // H2O2, LiqN2 Tanks and the Engine:                                     //
    //-----------------------------------------------------------------------//
    // Both Tanks are assume to be Circular Tori. The Major Diameter of both
    // Tanks is assumed to be "MaxD".
    // XXX: We currently assume that the LiqN2 and H2O2 Tanks are installed
    // just under the FuelTank, w/o any gaps (according to StarSem drawings):
    //
    constexpr static Len      LiqN2TankMinoR = 0.13_m;  // XXX: Rough approx
    constexpr static Len      H2O2TankMinoR  = 0.21_m;  // ditto

    //-----------------------------------------------------------------------//
    // Over-All Length and Angles of the StrapOn Booster Block:              //
    //-----------------------------------------------------------------------//
    constexpr static Len      H = TopConeH + MainTrCH + TailCylEnclH;
    static_assert(H.ApproxEquals(19.6_m));

    // The angle between the Booster's Xi axis and the LV's X axis. This is the
    // angle of the Main TrCone:
    constexpr static double   TanAlpha    =
      double((MaxD - TopConeD) / (2.0 * MainTrCH));

    constexpr static double   CosAlpha    = 1.0 / SqRt(1.0 + Sqr(TanAlpha));
    constexpr static double   SinAlpha    = TanAlpha * CosAlpha;

    // For info only: The angle of the TopCone:
    constexpr static double   TanAlphaTop = double(TopConeD / 2.0 / TopConeH);
    static_assert(TanAlpha  < TanAlphaTop);

    // The difference between the Alpha and AlphaTop angles should approximately
    // be equal to the angle of the MidTrCone of Stage2, due to the geometry of
    // StrapOn Boosters juction with Stage2:
    static_assert
      (Abs(ATan(TanAlphaTop) - ATan(TanAlpha) - ATan(S2::TanAlphaMid)) < 0.011);

    // The lowest point of the Booster's "TailCylEncl" (but w/o the extending
    // Engine Nozzles). NB: we must take into account the TailCylEncl radius.
    // The LowX co-ord closely matches that of Stage2:
    constexpr static Len      TailLowX    =
      TopX - H * CosAlpha - MaxD / 2.0 * SinAlpha;
    static_assert(Abs(TailLowX - S2::TailEnclLowX) < 0.005_m);

    // RD-107A Height:
    constexpr static Len      EngineH     = 2.578_m;

    // The lowest point of the Engine Nozzles (used for Thrust Momenta computa-
    // tions is the same as for Stage2). XXX: Similar to Stage2,  we currently
    // disregard the X-coord difference between the Verniers and the Main Engi-
    // ne:
    constexpr static Len      NozzlesLowX = S2::NozzlesLowX;
    static_assert(NozzlesLowX < TailLowX);
    constexpr static Len      NozzlesExtH = TailLowX - NozzlesLowX;

    //=======================================================================//
    // Masses:                                                               //
    //=======================================================================//
    // EmptyMass: XXX: StarSem says 3784 kg:
    constexpr static Mass     EmptyMass    = 3815.0_kg;

    // Mass of the RD-107A engine (part of EmptyMass):
    constexpr static Mass     EngMass      = 1090.0_kg;

    // Masses of Fuel (Naftil) and Oxidiser. As for Stage3, FuelMass includes
    // extra 0.2% for the antifreeze (2-EtoxyEthanol):
    constexpr static Mass     FuelMass     = 11458.0_kg * 1.002;
                                                         // StarSem: 11260 (T1)
    constexpr static Mass     OxidMass     = 27903.0_kg; // StarSem: 27900
    constexpr static Mass     H2O2Mass     = 1212.0_kg;
    constexpr static Mass     N2Mass       = 265.0_kg;
    constexpr static Mass     FullMass     = EmptyMass + FuelMass + OxidMass +
                                             H2O2Mass  + N2Mass;
    // N2Mass is split into a Liquid and Gaseous Phases, the initial masses are:
    constexpr static Mass     LiqN2Mass0   = 256.0_kg;
    constexpr static Mass     GasN2Mass0   = 9.0_kg;
    static_assert(LiqN2Mass0 + GasN2Mass0 == N2Mass);

    // UnSpendable Remnants of the Fuel and Oxidiser in Stage2 at the engine
    // cut-off time.   NB: The Remnants are about 1% of the corresp initial
    // masses. NB: They are Technically UnSpendable, so they do NOT include the
    // "guarantee margins"(???):
    constexpr static Mass     FuelRem      = 215.0_kg * 1.002;
    constexpr static Mass     OxidRem      = 451.0_kg;
    constexpr static Mass     H2O2Rem      = 125.0_kg;
    constexpr static Mass     LiqN2Rem     = 47.0_kg;

    //=======================================================================//
    // RD-107A (14D22) Engine Performance:                                   //
    //=======================================================================//
    // Isp of the whole Engine (SL/Vac, sec): quite consistent:
    // 263.1/320.0  (LPRE.DE),
    // 263.3/320.2  (EnergoMash, Wikipedia, red.is-telecom.ru),
    // 262  /319    (StarSem);
    // 269  /333.14 (Kravchenko; the vals seem to be too high),
    // so assume the higher vals:
    constexpr static Time     IspEngSL     = 263.3_sec;
    constexpr static Time     IspEngVac    = 320.2_sec;

    // Thrust of the whole Engine (SL/Vac, tf): also quite consistent:
    // 85.56/104    (Wikipedia,  Kravchenko)
    // 85.5 /104.14 (StarSem,    Zheleznyakov)
    // 85.6 /104.0  (LPRE.DE,    EnergoMash, red.is-telecom.ru, an.Wikipedia)
    // 79  / 96     (RosCosmos.ru/2115: w/o  the Verniers).
    // Again, we must have the MassRate constant across Vac/SL, so the closest
    // is 85.6/104.14. Thus:
    //
    constexpr static Force    ThrustEngSL    = 85600.0_kg * g0;
    constexpr static MassRate EngineMR       = ThrustEngSL / (IspEngSL * g0);
    // EngineMR     = 325.1 kg/sec
    constexpr static Force    ThrustEngVac   = IspEngVac  * g0 * EngineMR;
    // ThrustEngVac = 104.1 tf

    // Now the separate Isp and Thrust vals for the Main Engine  and the Verni-
    // ers. Unlike Stage2, here we do not have reliable Isp vals for the Verni-
    // ers, so will have to match the Thrust vals alone. This means that separ-
    // ate MassRates for the Main Engine and the Verniers will not be available,
    // but they are not needed anyway. And unlike Stage2, there are only 2 Ver-
    // niers in each Stage1 Booster StrapOn Block:
    constexpr static Force    ThrustMainSL   = 79000.0_kg * g0;
    constexpr static Force    ThrustMainVac  = 96000.0_kg * g0;
    constexpr static Force    ThrustVernSL2  = ThrustEngSL    - ThrustMainSL;
    constexpr static Force    ThrustVernVac2 = ThrustEngVac   - ThrustMainVac;
    constexpr static Force    ThrustVernSL1  = ThrustVernSL2  / 2.0;
    constexpr static Force    ThrustVernVac1 = ThrustVernVac2 / 2.0;

    // Separate Fuel and Oxid Rates are obtained using the Oxid/Fuel Ratio which
    // we derive from the over-all Fuel and Oxid spendable masses:
    constexpr static double OxidFuelRatio =
      double((OxidMass - OxidRem) / (FuelMass - FuelRem));  // ~2.44

    constexpr static double   FuelPart = 1.0  /   (1.0 + OxidFuelRatio);
    constexpr static double   OxidPart = OxidFuelRatio * FuelPart;

    constexpr static MassRate FuelMR   = EngineMR * FuelPart;
    constexpr static MassRate OxidMR   = EngineMR * OxidPart;

    //-----------------------------------------------------------------------//
    // RD-107A Ignition Sequence:                                            //
    //-----------------------------------------------------------------------//
    using FT = FlightTime;

    // Let t0=0 be the LiftOff ("Contact Separation") instant. Ignition occurs
    // at t0-15 sec (approx):
    // RD-107A thrust increases in stages  ("Preliminary", "1st Intermediate",
    // "2nd Intermediate", "Main"), where the "2nd Intermediate" occurs right
    // at "t0" and the "Main" (FullThrust) occurs at t0+6 sec:
    //
    constexpr static Time   IgnAdvance      = 15.0_sec;
    constexpr static FT     FullThrustTime  = FT(6.0_sec);  // Aka "t1"

    // Fuel and Oxidiser Mass @ t0=0. XXX: we assume that prior to t0, we run at
    // approx 25% on average (in reality, it changes over the time):
    constexpr static double IgnThrottlLevel = 0.25;
    constexpr static Mass   FuelMass0       =
      FuelMass  - FuelMR  * IgnAdvance * IgnThrottlLevel;
    constexpr static Mass   OxidMass0       =
      OxidMass  - OxidMR  * IgnAdvance * IgnThrottlLevel;

    // Fuel and Oxidiser Mass @ "FullThrustTime" ("t1", shortly after Lift-Off).
    // Prior to that, we assume a 75% throttling level, otherwise the thrust may
    // be insufficient for lift-off:
    constexpr static double LiftOffThrottlLevel = 0.75;
    constexpr static Mass   FuelMass1       =
      FuelMass0 -
      FuelMR  * (FullThrustTime - SC::LiftOffTime) * LiftOffThrottlLevel;

    constexpr static Mass   OxidMass1       =
      OxidMass0 -
      OxidMR  * (FullThrustTime - SC::LiftOffTime) * LiftOffThrottlLevel;

    //-----------------------------------------------------------------------//
    // RD-107A ShutDown Sequence:                                            //
    //-----------------------------------------------------------------------//
    // XXX: We assume that the Main Engine and the Verniers are throttled to
    // the 75% level:
    constexpr static Time   ThrottlAdvance        = 5.7_sec;
    constexpr static FT       ThrottlTime         =
      SC::Stage1CutOffTime -  ThrottlAdvance;
    constexpr static FT       CutOffTime          = SC::Stage1CutOffTime;

    // The Propellant Mass spent during the ShutDown sequence is:
    constexpr static double ShutDownThrottlLevel  = 0.75;
    constexpr static Mass   ShutDownSpentPropMass =
      EngineMR * ShutDownThrottlLevel * ThrottlAdvance;

    // Thus, the Propellant Mass left for the FullThrust Mode:
    constexpr static Mass     FullThrustPropMass  =
      FuelMass1 + OxidMass1 - ShutDownSpentPropMass - FuelRem - OxidRem;
    static_assert(IsPos(FullThrustPropMass));

    // Fuel and Oxid Masses @ "ThrottlTime":
    constexpr static Mass     FuelMassT =
      FuelMass1 - FuelMR   * (ThrottlTime - FullThrustTime);
    constexpr static Mass     OxidMassT =
      OxidMass1 - OxidMR   * (ThrottlTime - FullThrustTime);

    // Between ThrottlTime and CutOffTime, the following MassRates apply:
    constexpr static MassRate FuelMRT   = FuelMR * ShutDownThrottlLevel;
    constexpr static MassRate OxidMRT   = OxidMR * ShutDownThrottlLevel;

    //-----------------------------------------------------------------------//
    // RD-107A BurnTimes:                                                    //
    //-----------------------------------------------------------------------//
    // Then the Max Duration  at FullThrust is:
    constexpr static Time  MaxFullThrustDur =
      FullThrustPropMass / EngineMR;

    // Then the Max Time of RD-107A operation from the LiftOff (t0=0) (NOT from
    // Ignition which occurs before the LiftOff) to the is the following. Simi-
    // lar to Stage2, this is the End-Time, not Duration:
    constexpr static FT    MaxBurnTime      =
      FullThrustTime + MaxFullThrustDur + ThrottlAdvance;

    // "CutOffTime" must be less than (but close to) the "MaxBurnTime":
    static_assert(CutOffTime < MaxBurnTime);

    // H2O2:
    // NB: In addition, there is a MassRate due to H2O2 burning   to drive the
    // TurboPumps; however, this MassRate does not formally participate in the
    // Thrust. We assume that H2O2 flow begins ~10 sec before the LiftOff time
    // (before that, there is a passive flow of Fuel and Oxid) and lasts until
    // the end of the "MaxBurnTime". XXX: We also assume a constant H2O2 rate,
    // which is not exactly true before the Main (FullThrust) mode,   but the
    // error is insignificant in this case:
    //
    constexpr static Time     H2O2Adv = 10.0_sec;
    constexpr static MassRate H2O2MR  =
      (H2O2Mass - H2O2Rem) / (MaxBurnTime + H2O2Adv - SC::LiftOffTime);
    static_assert(IsPos(H2O2MR));

    // H2O2 Masses at LiftOff (t0=0) and at "FullThrustTime" (t1):
    constexpr static Mass     H2O2Mass0 =
      H2O2Mass  - H2O2Adv  *  H2O2MR;

    constexpr static Mass     H2O2Mass1 =
      H2O2Mass0 - (FullThrustTime - SC::LiftOffTime) * H2O2MR;

    // So  the total MassRate (at Full Thrust) is:
    constexpr static MassRate FullMR  = EngineMR + H2O2MR;

    // We assume that vaporisation of LiqN2 is also linear over time. however,
    // unlike H2O2, N2 is not exhaused, only re-distributed over the Tank vol-
    // umes becoming available:
    //
    constexpr static MassRate LiqN2MR =
      (LiqN2Mass0 - LiqN2Rem) / (MaxBurnTime - SC::LiftOffTime);
    static_assert(IsPos(LiqN2MR));

    // For testing: Minimal Mass of the Spent Stage2 (with all Remnants at their
    // minimal physical levels):
    constexpr static Mass     MinEndMass =
      EmptyMass + FuelRem + OxidRem + H2O2Rem + N2Mass;

    //-----------------------------------------------------------------------//
    // Fuel, Oxid, H2O2 and LiqN2 Masses @ "CutOffTime":                     //
    //-----------------------------------------------------------------------//
    constexpr static Mass     FuelMassC  =
      FuelMassT - FuelMRT  * (CutOffTime - ThrottlTime);
    constexpr static Mass     OxidMassC  =
      OxidMassT - OxidMRT  * (CutOffTime - ThrottlTime);
    constexpr static Mass     LiqN2MassC =
      LiqN2Mass0 - LiqN2MR * (CutOffTime - SC::LiftOffTime);
    constexpr static Mass     H2O2MassC  =
      H2O2Mass0  - H2O2MR  * (CutOffTime - SC::LiftOffTime);

   // Checks:
    static_assert
      (FuelMassC > FuelRem && OxidMassC > OxidRem && LiqN2MassC > LiqN2Rem &&
       H2O2MassC > H2O2Rem);

    static_assert((FuelMassC + OxidMassC).ApproxEquals
                  (FuelMassT + OxidMassT - ShutDownSpentPropMass));
  private:
    //=======================================================================//
    // "MechElement"s: "Proto"s with Yet-UnKnown Masses:                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The "TopTopConeProto" (Empty):                                        //
    //-----------------------------------------------------------------------//
    constexpr static TrC TopTopConeProto =
      TrC
      (
        TT::UnDef(),
        TopX,     TopY,  TopZ,
        CosAlpha, SinAlpha,
        CosPsi,   SinPsi,
        0.0_m,    TopTopConeD, TopTopConeH,
        Density(0.0)
      );

    //-----------------------------------------------------------------------//
    // "OxidTank*Proto"s:                                                    //
    //-----------------------------------------------------------------------//
    // Top: a SpherSegm:
    constexpr static SpS OxidTankTopProto =
      SpS
      (
        TT::UnDef(),
        true,                   // Facing Up
        TopTopConeProto.GetLow()[0],
        TopTopConeProto.GetLow()[1],
        TopTopConeProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        TopTopConeD,  OxidTankTopH,
        Propellants::LOxDens
      );

    // Up: a TrCone which is the low part of the "TopCone" (ie below the
    // "TopTopCone"):
    constexpr static TrC OxidTankUpProto =
      TrC
      (
        TT::UnDef(),
        TopTopConeProto.GetLow()[0],
        TopTopConeProto.GetLow()[1],
        TopTopConeProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        TopTopConeD,  OxidTankUpD, OxidTankUpH,
        Propellants::LOxDens
      );

    // Main (Mid): another TrCone:
    constexpr static TrC OxidTankMidProto =
      TrC
      (
        TT::UnDef(),
        OxidTankUpProto.GetLow()[0],
        OxidTankUpProto.GetLow()[1],
        OxidTankUpProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        OxidTankUpD,  OxidTankMidD,  OxidTankMidH,
        Propellants::LOxDens
      );

    // Bottom: a SpherSegm:
    constexpr static SpS OxidTankBtmProto =
      SpS
      (
        TT::UnDef(),
        false,                  // Facing Down
        OxidTankMidProto.GetLow()[0],
        OxidTankMidProto.GetLow()[1],
        OxidTankMidProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        OxidTankMidD, OxidTankBtmH,
        Propellants::LOxDens
      );

    //-----------------------------------------------------------------------//
    // "InterTankProto" (Empty):                                             //
    //-----------------------------------------------------------------------//
    constexpr static TrC InterTankProto =
      TrC
      (
        TT::UnDef(),
        OxidTankMidProto.GetLow()[0],
        OxidTankMidProto.GetLow()[1],
        OxidTankMidProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        InterTankUpD, InterTankLoD, InterTankH,
        Propellants::RG1Dens
      );

    //-----------------------------------------------------------------------//
    // "FuelTank*Proto"s:                                                    //
    //-----------------------------------------------------------------------//
    // Top: a SpherSegm:
    constexpr static SpS FuelTankTopProto =
      SpS
      (
        TT::UnDef(),
        true,                   // Facing Up
        InterTankProto.GetLow()[0],
        InterTankProto.GetLow()[1],
        InterTankProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        FuelTankTopD, FuelTankTopH,
        Propellants::RG1Dens
      );

    // Mid: a TrCone:
    constexpr static TrC FuelTankMidProto =
      TrC
      (
        TT::UnDef(),
        InterTankProto.GetLow()[0],
        InterTankProto.GetLow()[1],
        InterTankProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        FuelTankTopD, MaxD,  FuelTankMidH,
        Propellants::RG1Dens
      );

    // Bottom: a SpherSegm:
    constexpr static SpS FuelTankBtmProto =
      SpS
      (
        TT::UnDef(),
        false,                  // Facing Down
        FuelTankMidProto.GetLow()[0],
        FuelTankMidProto.GetLow()[1],
        FuelTankMidProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        MaxD,         FuelTankBtmH,
        Propellants::RG1Dens
      );

    //-----------------------------------------------------------------------//
    // Tail Cylindrical Enclosure:                                           //
    //-----------------------------------------------------------------------//
    // Attached to the Low of "FuelTankMidProto". The mass includes AeroFin:
    //
    constexpr static TrC TailCylEnclProto =
      TrC
      (
        TT::UnDef(),
        FuelTankMidProto.GetLow()[0],
        FuelTankMidProto.GetLow()[1],
        FuelTankMidProto.GetLow()[2],
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        MaxD,         MaxD,     TailCylEnclH,
        Density(0.0)            // Does not contain any Propellants by itself
      );

    //-----------------------------------------------------------------------//
    // "LiqN2Tank*" Protos:                                                  //
    //-----------------------------------------------------------------------//
    // NB: Unlike Stage2, in a Stage1 Booster the LiqN2Tank is located ABOVE the
    // H2O2Tank, and has a slightly smaller diameter than the latter:
    //
    constexpr static Len LiqN2TankBaseCentrX =
        FuelTankBtmProto.GetLow()[0] - LiqN2TankMinoR * CosAlpha;

    constexpr static Len LiqN2TankBaseCentrY =
        FuelTankBtmProto.GetLow()[1] + LiqN2TankMinoR * SinAlpha * CosPsi;

    constexpr static Len LiqN2TankBaseCentrZ =
        FuelTankBtmProto.GetLow()[2] + LiqN2TankMinoR * SinAlpha * SinPsi;

    constexpr static Len LiqN2TankD          = MaxD - 2.0 * H2O2TankMinoR;

    constexpr static Tor LiqN2TankTopProto   =
      Tor
      (
        TT::UnDef(),
        true,                   // Facing Up
        LiqN2TankBaseCentrX, LiqN2TankBaseCentrY, LiqN2TankBaseCentrZ,
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        LiqN2TankMinoR * 2.0,
        LiqN2TankMinoR,
        LiqN2TankD,
        Propellants::LN2Dens
      );

    constexpr static Tor LiqN2TankBtmProto =
      Tor
      (
        TT::UnDef(),
        false,                  // Facing Down
        LiqN2TankBaseCentrX, LiqN2TankBaseCentrY, LiqN2TankBaseCentrZ,
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        LiqN2TankMinoR * 2.0,
        LiqN2TankMinoR,
        LiqN2TankD,
        Propellants::LN2Dens
      );

    //-----------------------------------------------------------------------//
    // "H2O2Tank*"  Protos:                                                  //
    //-----------------------------------------------------------------------//
    constexpr static Len H2O2TankBaseCentrX =
        LiqN2TankBtmProto.GetLow()[0] - H2O2TankMinoR * CosAlpha;

    constexpr static Len H2O2TankBaseCentrY =
        LiqN2TankBtmProto.GetLow()[1] + H2O2TankMinoR * SinAlpha * CosPsi;

    constexpr static Len H2O2TankBaseCentrZ =
        LiqN2TankBtmProto.GetLow()[2] + H2O2TankMinoR * SinAlpha * SinPsi;

    constexpr static Tor H2O2TankTopProto  =
      Tor
      (
        TT::UnDef(),
        true,                   // Facing Up
        H2O2TankBaseCentrX, H2O2TankBaseCentrY, H2O2TankBaseCentrZ,
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        H2O2TankMinoR * 2.0,
        H2O2TankMinoR,
        MaxD,
        Propellants::H2O2Dens
      );

    constexpr static Tor H2O2TankBtmProto =
      Tor
      (
        TT::UnDef(),
        false,                  // Facing Down
        H2O2TankBaseCentrX, H2O2TankBaseCentrY, H2O2TankBaseCentrZ,
        CosAlpha,     SinAlpha,
        CosPsi,       SinPsi,
        H2O2TankMinoR * 2.0,
        H2O2TankMinoR,
        MaxD,
        Propellants::H2O2Dens
      );

    //-----------------------------------------------------------------------//
    // Scale Factor to determine the individual masses:                      //
    //-----------------------------------------------------------------------//
    // XXX: Once again, InterStageGrid is not modeled, its mass is spread over
    // the masses of other "ME"s:
    //
    constexpr static double ScaleFactor =
      ME::GetMassScale
      (
        // Components for which we provide the TotalMass (below):
        { &TopTopConeProto,
          &OxidTankTopProto,  &OxidTankUpProto,  &OxidTankMidProto,
          &OxidTankBtmProto,
          &InterTankProto,
          &FuelTankTopProto,  &FuelTankMidProto, &FuelTankBtmProto,
          &TailCylEnclProto,
          &LiqN2TankTopProto, &LiqN2TankBtmProto,
          &H2O2TankTopProto,  &H2O2TankBtmProto
        },
        EmptyMass - EngMass
      );

  public:
    //=======================================================================//
    // "MechElement"s with Real Masses:                                      //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The "TopTopCone": Empty:                                              //
    //-----------------------------------------------------------------------//
    constexpr static TrC TopTopCone   =
      ME::ProRateMass(TopTopConeProto,   ScaleFactor);

    //-----------------------------------------------------------------------//
    // OxidTank-Related:                                                     //
    //-----------------------------------------------------------------------//
    constexpr static SpS OxidTankTop  =
      ME::ProRateMass(OxidTankTopProto,  ScaleFactor);

    constexpr static TrC OxidTankUp   =
      ME::ProRateMass(OxidTankUpProto,   ScaleFactor);

    constexpr static TrC OxidTankMid  =
      ME::ProRateMass(OxidTankMidProto,  ScaleFactor);

    constexpr static SpS OxidTankBtm  =
      ME::ProRateMass(OxidTankBtmProto,  ScaleFactor);

    //-----------------------------------------------------------------------//
    // InterTank:                                                            //
    //-----------------------------------------------------------------------//
    constexpr static TrC InterTank    =
      ME::ProRateMass(InterTankProto,    ScaleFactor);

    // There must be a Gap between OxidTankBtm and FuelTankTop:
    constexpr static Len OxidFuelTankGap =
      (OxidTankBtmProto.GetLow()[0] - FuelTankTopProto.GetUp()[0]) / CosAlpha;
    static_assert(IsPos(OxidFuelTankGap));

    //-----------------------------------------------------------------------//
    // FuelTank-Related:                                                     //
    //-----------------------------------------------------------------------//
    constexpr static SpS FuelTankTop  =
      ME::ProRateMass(FuelTankTopProto,  ScaleFactor);

    constexpr static TrC FuelTankMid  =
      ME::ProRateMass(FuelTankMidProto,  ScaleFactor);

    constexpr static SpS FuelTankBtm  =
      ME::ProRateMass(FuelTankBtmProto,  ScaleFactor);

    //-----------------------------------------------------------------------//
    // Tail Section:                                                         //
    //-----------------------------------------------------------------------//
    // Tail Cylindrical Enclosure:
    constexpr static TrC TailCylEncl  =
      ME::ProRateMass(TailCylEnclProto,  ScaleFactor);

    // LiqN2Tank-Related:
    constexpr static Tor LiqN2TankTop =
      ME::ProRateMass(LiqN2TankTopProto, ScaleFactor);

    constexpr static Tor LiqN2TankBtm =
      ME::ProRateMass(LiqN2TankBtmProto, ScaleFactor);

    // H2O2Tank-Related:
    constexpr static Tor H2O2TankTop =
      ME::ProRateMass(H2O2TankTopProto,  ScaleFactor);

    constexpr static Tor H2O2TankBtm =
      ME::ProRateMass(H2O2TankBtmProto,  ScaleFactor);

    // Verify the lowest X-coord on the axis of the "TailCylEncl":
    static_assert(TailCylEncl.GetLow()[0].ApproxEquals(TopX - H * CosAlpha));

    // FIXME: It looks like, with the curr geometry, there is a sufficiently
    // large usused space between the bottom of the H2O2Tank and the Engine.
    // It is not clear whether it is real, or it is a result of our imprecise
    // geometry. It is also not clear whether "EngineH" includes the support
    // structure...
    constexpr static Len TailGapH =
      H2O2TankBtm.GetLow()[0] - EngineH - NozzlesLowX;
    static_assert(IsPos(TailGapH));       // ~1.02 m!

    // RD-107A Engine, modeled as a PointMass.   XXX: We assume that its CoM is
    // located at 1/3 height from the Engine Top (ie 2/3 height from Nozzles),
    // but this is obviously GROSSLY IMPRECISE.
    // NB: the Engine is installed parallel to Stage2 Engine, NOT parallel to
    // the Booster OX axis, so the Alpha angle does not apply in this case:
    //
    constexpr static Len EngineCoMX = NozzlesLowX + (2.0/3.0) * EngineH;
    constexpr static Len EngineCoMY = H2O2TankBtm.GetLow()[1];
    constexpr static Len EngineCoMZ = H2O2TankBtm.GetLow()[2];

    constexpr static PM  Engine     =
      PM
      (
        TT::UnDef(),
        EngineCoMX, EngineCoMY, EngineCoMZ,
        EngMass
      );

    //-----------------------------------------------------------------------//
    // Empty Stage1 Booster as a whole:                                      //
    //-----------------------------------------------------------------------//
    // Similar to Stage2 and dissimilar to Stage3, the pressurisation gas (N2),
    // although not expendable, NOT included into the "EmptyME":
    //
    constexpr static ME  EmptyME  =
      TopTopCone   +
      OxidTankTop  + OxidTankUp   + OxidTankMid + OxidTankBtm +
      InterTank    +
      FuelTankTop  + FuelTankMid  + FuelTankBtm +
      TailCylEncl  +
      LiqN2TankTop + LiqN2TankBtm +
      H2O2TankTop  + H2O2TankBtm  +
      Engine;

    // The over-all EmptyMass check:
    static_assert(EmptyME.GetMass().ApproxEquals(EmptyMass));

  private:
    //=======================================================================//
    // Propellant Volumes and Mass Capacities:                               //
    //=======================================================================//
    // Oxid:
    constexpr static Mass OxidTankTopMC      = OxidTankTop.GetPropMassCap();
    constexpr static Mass OxidTankUpMC       = OxidTankUp .GetPropMassCap();
    constexpr static Mass OxidTankMidMC      = OxidTankMid.GetPropMassCap();
    constexpr static Mass OxidTankBtmMC      = OxidTankBtm.GetPropMassCap();
    // Unions:
    constexpr static Mass OxidTankBtmMidMC   = OxidTankBtmMC    + OxidTankMidMC;
    constexpr static Mass OxidTankBtmMidUpMC = OxidTankBtmMidMC + OxidTankUpMC;

    // Fuel:
    constexpr static Mass FuelTankTopMC      = FuelTankTop.GetPropMassCap();
    constexpr static Mass FuelTankMidMC      = FuelTankMid.GetPropMassCap();
    constexpr static Mass FuelTankBtmMC      = FuelTankBtm.GetPropMassCap();
    // Union:
    constexpr static Mass FuelTankBtmMidMC   = FuelTankBtmMC    + FuelTankMidMC;

    // LiqN2:
    constexpr static Mass LiqN2TankTopMC     = LiqN2TankTop.GetPropMassCap();
    constexpr static Mass LiqN2TankBtmMC     = LiqN2TankBtm.GetPropMassCap();

    // H2O2:
    constexpr static Mass H2O2TankTopMC      = H2O2TankTop.GetPropMassCap();
    constexpr static Mass H2O2TankBtmMC      = H2O2TankBtm.GetPropMassCap();

  public:
    // Maximum Theoretical Capacities of the resp Tanks:
    constexpr static Mass OxidTankMC   = OxidTankBtmMidUpMC + OxidTankTopMC;
    constexpr static Mass FuelTankMC   = FuelTankBtmMidMC   + FuelTankTopMC;
    constexpr static Mass LiqN2TankMC  = LiqN2TankTopMC     + LiqN2TankBtmMC;
    constexpr static Mass H2O2TankMC   = H2O2TankTopMC      + H2O2TankBtmMC;

    // Propellant Load Ratios (ActualMass  / TheorMassCapacity):
    constexpr static double OxidLoadRatio  = double(OxidMass   / OxidTankMC);
    static_assert(OxidLoadRatio  < 1.0);
    constexpr static double FuelLoadRatio  = double(FuelMass   / FuelTankMC);
    static_assert(FuelLoadRatio  < 1.0);
    constexpr static double LiqN2LoadRatio = double(LiqN2Mass0 / LiqN2TankMC);
    static_assert(LiqN2LoadRatio < 1.0);
    constexpr static double H2O2LoadRatio  = double(H2O2Mass   / H2O2TankMC);
    static_assert(H2O2LoadRatio  < 1.0);

  private:
    //-----------------------------------------------------------------------//
    // "ME" objs for Max Theoretical Propellant Loads in Tank Sections:      //
    //-----------------------------------------------------------------------//
    // (For optimisation of "GetDynParams"):
    // Oxid:
    constexpr static ME OxidTopME         = OxidTankTop.GetPropBulkME();
    constexpr static ME OxidUpME          = OxidTankUp .GetPropBulkME();
    constexpr static ME OxidMidME         = OxidTankMid.GetPropBulkME();
    constexpr static ME OxidBtmME         = OxidTankBtm.GetPropBulkME();
    // Unions:
    constexpr static ME OxidBtmMidME      = OxidBtmME      + OxidMidME;
    constexpr static ME OxidBtmMidUpME    = OxidBtmMidME   + OxidUpME;
    // The whole maximum Oxid bulk:
    constexpr static ME OxidME            = OxidBtmMidUpME + OxidTopME;

    //Fuel:
    constexpr static ME FuelTopME         = FuelTankTop.GetPropBulkME();
    constexpr static ME FuelMidME         = FuelTankMid.GetPropBulkME();
    constexpr static ME FuelBtmME         = FuelTankBtm.GetPropBulkME();
    // Union:
    constexpr static ME FuelBtmMidME      = FuelBtmME      + FuelMidME;
    // The whole maximum Fuel bulk:
    constexpr static ME FuelME            = FuelBtmMidME   + FuelTopME;

    // For LiqN2 and H2O2, "Top" MEs are not required:
    // LiqN2:
    constexpr static ME LiqN2BtmME        = LiqN2TankBtm.GetPropBulkME();
    // H2O2:
    constexpr static ME H2O2BtmME         = H2O2TankBtm .GetPropBulkME();

  public:
    //=======================================================================//
    // Flight Control:                                                       //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Thrust Vector Control:                                                //
    //-----------------------------------------------------------------------//
    // Achieved via 2 Gimbaled Vernier Chambers. A Pair of Vernier Chambers is
    // located (notionally) in the XY or XZ plane. They are gimbaled SIMULTANE-
    // OUSLY (NOT independently) within the GimbalAmpl in the Tangential Plane.
    // Similar to Stages 3 and 2, we assume that Gimbaling (Deflection) Angle is
    // positive when the corresp Verniers are deflected Counter-Clock-Wise wrt
    // the corresp rotation axis (same as the "attachment direction" of this
    // Block):
    //   B: @ -Y (Angle > 0 ==> chamber towards -Z, thrust towards +Z);
    //   V: @ -Z (Angle > 0 ==> chanber towards +Y, thrust towards -Y);
    //   G: @ +Y (Angle > 0 ==> chamber towards +Z, thrust towards -Z);
    //   D: @ +Z (Angle > 0 ==> chamber towards -Y, thrust towards +Y).
    //
    // When viewed from +X, positive deflection angle corredponds to a Counter-
    // Clock-Wise movement of the Nozzles projections in the YZ plane.   Thus,
    // gimbaling of a Booster Block is characterised by just 1 angle with the
    // following amplitude:
    //
    constexpr static Angle_deg VernGimbalAmpl = Angle_deg(45.0);

    //-----------------------------------------------------------------------//
    // AeroFin Control:                                                      //
    //-----------------------------------------------------------------------//
    // Each AeroFin is approx a right equilateral triangle with Sides = 1 m:
    constexpr static Area      AeroFinS       = Area(0.5);

    // XXX: The max AeroFin rotation angle is unknown, assume 30 deg (from an
    // available photo).  Larger values would probably cause too large angles
    // of attack.
    // AeroFin is rotatable along the same "block attachment direction" as used
    // for Verniers Gimbaling (that is, -Y, -Z, +Y or +Z).  Rotation angle is
    // positive when rotation is performed Counter-Clock-Wise as viewed  from
    // the corresp end of the rotation axis (similar to Verniers deflection):
    //
    constexpr static Angle_deg AeroFinAmpl    = Angle_deg(30.0);

    //=======================================================================//
    // Dynamic Params as a function of Flight Time:                          //
    //=======================================================================//
    // NB: This method is NOT "constexpr": it is intended to be called at Run-
    // Time (eg multiple times during Trajectory Integration).  Because Stage1
    // Booster has an AeroFin control in addition to the Verniers, we need the
    // velocity vector to be provided as an arg:
    //
    static StageDynParams<LVSC::Soyuz21b>
    GetDynParams
    (
      FT               a_ft,
      Pressure         a_p,           // Curr Atmospheric Pressure
      ME::VelVE const& a_v,           // In the ECOS
      Angle_deg        a_verns_defl,  // Same angle for both Verns
      Angle_deg        a_aerofin_defl // AeroFin deflection angle
    );
  };

  //=========================================================================//
  // Stage1 Blocks ('B', 'V', 'G', 'D'):                                     //
  //=========================================================================//
  using Soyuz21b_BlockB = Soyuz21b_Stage1_Booster<'B'>;
  using Soyuz21b_BlockV = Soyuz21b_Stage1_Booster<'V'>;
  using Soyuz21b_BlockG = Soyuz21b_Stage1_Booster<'G'>;
  using Soyuz21b_BlockD = Soyuz21b_Stage1_Booster<'D'>;
}
// End namespace SpaceBallistivs
