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

    // The X-coord of the Booster top: Relative to MaxD of Stage2:
    static_assert(S2::OxidTankUp.GetLow()[0] == S2::OxidTankLow.GetUp()[0]);
    constexpr static Len      TopX   = S2::OxidTankUp.GetLow()[0] - 0.56_m;

    // XXX: The following geometry of is somewhat approximate, it does not take
    // into account the full complicated shape of the Booster Block...

    // The TopCone (assumed to be a full cone, not truncated!):
    constexpr static Len      TopConeH     = 4.015_m;
    constexpr static Len      TopConeD     = 1.350_m;

    // The Top-most part of the TopCone is by itself empty, and is capped with a
    // ball joint used for attaching the StrapOn Booster to the core Stage2:
    constexpr static double   TopTopConeF  = 0.38;
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
    constexpr static Len      MainTrCH     = 11.165_m;
    constexpr static Len      MaxD         = 2.680_m;
    constexpr static Len      OxidTankMidH = 8.0_m;
    constexpr static Len      OxidTankMidD =
      OxidTankUpD  + (MaxD - OxidTankUpD)  * double(OxidTankMidH / MainTrCH);

    // Below is the InterTank section which is empty by itself and serves as an
    // enclosure for the OxidTankBtm and FuelTankTop SpherSegms:
    constexpr static Len      InterTankH   = 1.075_m;
    constexpr static Len      InterTankUpD = OxidTankMidD;
    constexpr static Len      InterTankLoD =
      InterTankUpD + (MaxD - OxidTankUpD)  * double(InterTankH   / MainTrCH);

    // OxidTankBtm:
    constexpr static Len      OxidTankBtmH = 0.45_m;
    constexpr static Len      OxidTankBtmD = InterTankUpD;

    //-----------------------------------------------------------------------//
    // FuelTank:                                                             //
    //-----------------------------------------------------------------------//
    // FuelTankTop:
    constexpr static Len      FuelTankTopH = 0.48_m;
    static_assert(IsPos(InterTankH - OxidTankBtmH - FuelTankTopH));
    constexpr static Len      FuelTankTopD = InterTankLoD;

    // FuelTankMid: Its bottom diameter is "MaxD":
    constexpr static Len      FuelTankMidH =
      MainTrCH - OxidTankMidH - InterTankH;
    static_assert(IsPos(FuelTankMidH));

    // Then the cylindrical enclosure of the FuelTankBtm, H2O2Tank, LiqN2 Tank
    // and the Engine, of the "MaxD" diameter. XXX: We ignore the fact that it
    // is not strictly cylindrical in its lowest (Engine) part:
    //
    constexpr static Len      TailCylH     = 4.020_m;
    // FuelTankBtm:
    constexpr static Len      FuelTankBtmH = 0.52_m;

    //-----------------------------------------------------------------------//
    // H2O2, LiqN2 Tanks and the Engine:                                     //
    //-----------------------------------------------------------------------//
    // TODO...

    //-----------------------------------------------------------------------//
    // Over-All Length and Angles of the StrapOn Booster Block:              //
    //-----------------------------------------------------------------------//
    constexpr static Len      H = TopConeH + MainTrCH + TailCylH;
    static_assert(H.ApproxEquals(19.2_m));

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
      (Abs(ATan(TanAlphaTop) - ATan(TanAlpha) - ATan(S2::TanAlphaMid)) < 0.01);

    // The lowest point of the Booster's "TailCyl" (but w/o the extending Engine
    // Nozzles):
    constexpr static Len      TailLowX    = TopX - H * CosAlpha;

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
    // Let t0=0 be the LiftOff ("Contact Separation") instant. Ignition occurs
    // at t0-15 sec (approx):
    // RD-107A thrust increases in stages  ("Preliminary", "1st Intermediate",
    // "2nd Intermediate", "Main"), where the "2nd Intermediate" occurs right
    // at "t0" and the "Main" (FullThrust) occurs at t0+6 sec:
    //
    constexpr static Time   IgnAdvance      = 15.0_sec;
    constexpr static Time   FullThrustTime  =  6.0_sec;   // Aka "t1"

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
      FuelMass0 - FuelMR  * FullThrustTime  * LiftOffThrottlLevel;
    constexpr static Mass   OxidMass1       =
      OxidMass0 - OxidMR  * FullThrustTime  * LiftOffThrottlLevel;

    //-----------------------------------------------------------------------//
    // RD-107A ShutDown Sequence:                                            //
    //-----------------------------------------------------------------------//
    // XXX: We assume that the Main Engine and the Verniers are throttled to
    // the 75% level:
    constexpr static Time   ThrottlAdvance        = 5.7_sec;
    constexpr static Time     ThrottlTime         =
      SC::Stage1CutOffTime -  ThrottlAdvance;
    constexpr static Time     CutOffTime          = SC::Stage1CutOffTime;

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
    constexpr static Time  MaxBurnTime      =
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
      (H2O2Mass - H2O2Rem) / (H2O2Adv + MaxBurnTime);
    static_assert(IsPos(H2O2MR));

    // H2O2 Masses at LiftOff (t0=0) and at "FullThrustTime" (t1):
    constexpr static Mass     H2O2Mass0 = H2O2Mass  - H2O2Adv        * H2O2MR;
    constexpr static Mass     H2O2Mass1 = H2O2Mass0 - FullThrustTime * H2O2MR;

    // So  the total MassRate (at Full Thrust) is:
    constexpr static MassRate FullMR  = EngineMR + H2O2MR;

    // We assume that vaporisation of LiqN2 is also linear over time. however,
    // unlike H2O2, N2 is not exhaused, only re-distributed over the Tank vol-
    // umes becoming available:
    //
    constexpr static MassRate LiqN2MR = (LiqN2Mass0 - LiqN2Rem) / MaxBurnTime;
    static_assert(IsPos(LiqN2MR));

    // For testing: Minimal Mass of the Spent Stage2 (with all Remnants at their
    // minimal physical levels):
    constexpr static Mass     MinEndMass =
      EmptyMass + FuelRem + OxidRem + H2O2Rem + N2Mass;

    //-----------------------------------------------------------------------//
    // Fuel, Oxid, H2O2 and LiqN2 Masses @ "CutOffTime":                     //
    //-----------------------------------------------------------------------//
    constexpr static Mass     FuelMassC =
      FuelMassT - FuelMRT * (CutOffTime - ThrottlTime);
    constexpr static Mass     OxidMassC =
      OxidMassT - OxidMRT * (CutOffTime - ThrottlTime);
    constexpr static Mass    H2O2MassC  = H2O2Mass0  - H2O2MR  * CutOffTime;
    constexpr static Mass    LiqN2MassC = LiqN2Mass0 - LiqN2MR * CutOffTime;

   // Checks:
    static_assert
      (FuelMassC  > FuelRem && OxidMassC > OxidRem && H2O2MassC > H2O2Rem &&
       LiqN2MassC > LiqN2Rem);

    static_assert((FuelMassC + OxidMassC).ApproxEquals
                  (FuelMassT + OxidMassT - ShutDownSpentPropMass));

  public:
    //=======================================================================//
    // Flight Control:                                                       //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Thrust Vector Control:                                                //
    //-----------------------------------------------------------------------//
    // Achieved via 2 Gimbaled Vernier Chambers. A Pair of Vernier Chambers is
    // located (notionally) in the XY or XZ plane:
    //   B: @ -Y (gimbaled in Z)
    //   V: @ -Z (gimbaled in Y)
    //   G: @ +Y (gimbaled in Z)
    //   D: @ +Z (gimbaled in Y)
    // Each Vernier can be deflected within the GimbalAmpl in the Tangential
    // Plane. Similar to Stage3 and Stage2, we assume that Gimbaling (Deflect-
    // ion) Angle s positive when the corresp Vernier is deflected  Counter-
    // Clock-Wise in the YZ plane, as viewed from the positive end of the OX
    // axis. And again, similar to Stages 3 and 2,   opposite Verniers would
    // normally be moved symmetrically in one direction,
    // so that
    //    VernDeflections[0] = - VernDeflections[2],
    //    VernDeflections[1] = - VernDeflections[3],
    // but this is not strictly enforced:
    //
    constexpr static Angle_deg VernGimbalAmpl = Angle_deg(45.0);
    using            VernDeflections          = std::array<Angle_deg, 4>;

    //-----------------------------------------------------------------------//
    // AeroFin Control:                                                      //
    //-----------------------------------------------------------------------//
    // Each AeroFin is approx a right equilateral triangle with Sides = 1 m:
    constexpr static Area      AeroFinS       = Area(0.5);

    // XXX: The max AeroFin rotation angle is unknown, assume 30 deg (from an
    // available photo).  Larger values would probably cause too large angles
    // of attack:
    constexpr static Angle_deg AeroFinAmpl    = Angle_deg(30.0);

    //=======================================================================//
    // Dynamic Params as a function of Flight Time:                          //
    //=======================================================================//
    // NB: This method is NOT "constexpr": it is intended to be called at Run-
    // Time (eg multiple times during Trajectory Integration):
    //
    static StageDynParams<LVSC::Soyuz21b>
    GetDynParams
    (
      Time                   a_t,
      Pressure               a_p,      // Curr Atmospheric Pressure
      VernDeflections const& a_vern_defls,
      Angle_deg              a_aerofin_defl
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
