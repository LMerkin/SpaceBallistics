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
    // Block B: +Y, Psi=0
    // Block V: +Z, Psi=Pi/2
    // Block G: -Y, Psi=Pi
    // Block D: -Z, Psi-3*Pi/2
    //
    // The Psi angle in the OYZ plane (see "RotationShell" for the definition):
    // Psi = Pi/2 * I:
    constexpr static double   CosPsi =
      (Block == 'B') ? 1.0 : (Block == 'G') ? -1.0 : 0.0;
    constexpr static double   SinPsi =
      (Block == 'V') ? 1.0 : (Block == 'D') ? -1.0 : 0.0;

    // The X-coord of the Booster top: Relative to MaxD of Stage2:
    static_assert(S2::OxidTankUp.GetLow()[0] == S2::OxidTankLow.GetUp()[0]);
    constexpr static Len      TopX   = S2::OxidTankUp.GetLow()[0] - 0.56_m;

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
    // XXX: For Vernier Chambers, use the same vals as in Stage2, but unlike
    // Stage2, here there are only 2 Verniers per block,  deflectable in the
    // Tangential plane:
    constexpr static Time     IspVernSL1     = 251.9_sec;  // LPRE.DE
    constexpr static Time     IspVernVac1    = 313.1_sec;  //
    constexpr static Force    ThrustVernSL1  = 2700.0_kg * g0;
    constexpr static Force    ThrustVernVac1 =             // 3.356 tf
      ThrustVernSL1 * double(IspVernVac1 / IspVernSL1);
    constexpr static Force    ThrustVernSL2  = 2.0 * ThrustVernSL1;
    constexpr static Force    ThrustVernVac2 = 2.0 * ThrustVernVac1;

    // Isp of the Main Engine (SL/Vac, sec):
    // 263.1/320.0 (LPRE.DE), 263.3/320.2 (EnergoMash etc), 262/319 (StarSem);
    // similar to Stage2, assume the higher values for the Main Engine   (and
    // slightly lower vals for the Vernier Chambers):
    //
    constexpr static Time     IspMainSL      = 263.3_sec;
    constexpr static Time     IspMainVac     = 320.2_sec;

    // Thrust of the Main Engine (SL/Vac, tf):
    // 85.5/104.14 (StarSem, perhaps incl the Verniers?),
    // 79.6/ 97.0  (StarSem,         excl the Verniers -- computed),
    // 85.6/104.0  (EnergoMash,      incl the Verniers?)
    //
    // We will use the EnergoMash SL value, adjust it for the Verniers, and
    // pro-rate to get the exact Vac value:
    //
    constexpr static Force    ThrustMainSL   =
      85600.0_kg * g0 - ThrustVernSL2;                      // 80.2  tf
    constexpr static Force    ThrustMainVac  =
      ThrustMainSL * double(IspMainVac / IspMainSL);        // 97.53 tf

    // We can then calculate the MassRates for the Main Engine and for each
    // Vernier Chamber:
    constexpr static MassRate MainMR   = ThrustMainSL  / (IspMainSL  * g0);
    static_assert(MainMR.ApproxEquals   (ThrustMainVac / (IspMainVac * g0)));
                                                            // 304.6 kg/sec
    constexpr static MassRate VernMR1  = ThrustVernSL1 / (IspVernSL1 * g0);
    static_assert
      (VernMR1.ApproxEquals(ThrustVernVac1 / (IspVernVac1 * g0)));
    // ~10.72  kg/sec, NOT 4.15+8.55=12.70 kg/sec as it would be with the orig-
    // inal RD-107 Vernier Chamber data...

    constexpr static MassRate VernMR2  = 2.0 * VernMR1;

    // Total MassRate for the whole Engine at FullThrust:
    constexpr static MassRate EngineMR = MainMR + VernMR2;  // ~326.04 kg/sec

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
