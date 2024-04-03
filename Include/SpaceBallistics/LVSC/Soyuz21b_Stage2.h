// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/LVSC/Soyuz21b_Stage2.h":                 //
//         Mathematical Model of the "Soyuz-2.1b" Stage2 ("Block A")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/MechElement.hpp"
#include "SpaceBallistics/Soyuz21b_Consts.h"
#include "SpaceBallistics/Soyuz21b_Head.h"
#include "SpaceBallistics/Propellants.h"
#include "SpaceBallistics/StageDynParams.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC = Soyuz21b_Consts;
  namespace SH = Soyuz21b_Head;
  using     ME = MechElement;

  //=========================================================================//
  // "Soyuz21b_Stage2" Class:                                                //
  //=========================================================================//
  // The following local co-ord system is used for convenience (so that most
  // X-coords are positive and the origin  is not  affected by separation of
  // stages):
  // (*) The OX axis is the main axis of the rocket. The OX positive direction
  //     is towards the TAIL (REAL of the LV). The origin O is the LOWER base
  //     of the InterStage (junction plane with Stage2).  That is, X <= 0 for
  //     for InterStage, the optional Stage4, Payload Adapter/Dispenser, Pay-
  //     load itself and the Fairing, and X >= 0 for Stages 3, 2, 1.
  // (*) The OY and OZ axes are such that the OXY and OXZ planes pass through
  //     the symmetry axes of the corresp opposite strap-on boosters (Blocks
  //     B, V, G, D -- Stage 1), and OXYZ is a right-oriented co-ords system.
  // (*) This class is actually a namespace: all entities declared in it are
  //     "static" (and mostly "constexpr") ones. However, we still make it a
  //     class, not a namespace, in order to be able to control the external
  //     visibility of its entities:
  //
  class Soyuz21b_Stage2
  {
  private:
    //=======================================================================//
    // No objects can be created of this class:                              //
    //=======================================================================//
    Soyuz21b_Stage2() = delete;

  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Geometry:                                                             //
    //-----------------------------------------------------------------------//
    // Stage2 Top (Base of Deflector Cone): Bottom of Stage3 minus the Height
    // of InterStage Grid:
    constexpr static Len    GridH         = 0.975_m;
    constexpr static Len    TopX          = SC::X0 - SC::Stage3Len - GridH;

    // Deflector Code Height and Base Diameter:
    constexpr static Len    DeflConeH     = GridH / 2.0;
    constexpr static Len    TopD          = 2.700_m;  // Slightly >than Stage3D

    // Equipment Bay: Cylindrical or expanding TrCone? Information differs in
    // different drawings. Assume Cylindrical, of the same "DeflConeD" diameter:
    constexpr static Len    EquipBayH     = 0.595_m;

    // The over-all Upper (Expanding) TrCone: Contains the Top SpherSegm of the
    // OxidTank, and the Upper TrCone of the OxidTank itself. The diameter inc-
    // reases from "DeflConeD" to "MaxD":
    constexpr static Len    UpperTrCH     = 5.835_m;
    constexpr static Len    MaxD          = 2.950_m;  // Max for entire Stage2

    // The over-all Middle (Narrowing-Down) TrCone: Contains the Low TrCone and
    // Bottom SpherSegm of the OxidTank, as well as Top SpherSegm of the Fuel
    // Tank. The diameter decreses from MaxD to MinD:
    constexpr static Len    MidTrCH       = 4.55_m;
    constexpr static Len    MinD          = 2.05_m;

    // The Lower Cylinder: Contains the Mid Cylinder and the Botrtom ShperSegm
    // of the FuelTank, the H2O2 Tank and the N2 (LiqNit) Tank (see below).

    // OxidTank:
    // Top SpherSegm:
    constexpr static Len    OxidTankTopH  = 0.900_m;  // Approx
    constexpr static Len    OxidTankTopD  =
      TopD + double(OxidTankTopH / UpperTrCH) * (MaxD - TopD);

    // Height of the Upper TrCone part of the OxidTank:
    constexpr static Len    OxidTankUpH   = UpperTrCH - OxidTankTopH;

    // Height of the Lower TrCone part of the OxidTank:
    constexpr static Len    OxidTankLowH  = 3.20_m;

    // Height and Base Diameter of the Bottom SpherSegm part of the OxidTank:
    constexpr static Len    OxidTankBtmH  = 0.635_m;
    constexpr static Len    OxidTankBtmD  =
      MaxD + double(OxidTankLowH / MidTrCH) * (MinD - MaxD);

    // FuelTank:
    // Top SpherSegm (Base Diameter = MinD):
    constexpr static Len    FuelTankTopH  = 0.5_m;

    // There is a gap between the OxidTankBtm and FuelTankTop:
    constexpr static Len    InterBankBayH =
      MidTrCH - OxidTankLowH - OxidTankBtmH - FuelTankTopH;
    static_assert(IsPos(InterBankBayH));

    // Length of the main cylindrical section of the FuelTank:
    constexpr static Len    FuelTankMidH  =

    //-----------------------------------------------------------------------//
    // Masses:                                                               //
    //-----------------------------------------------------------------------//
    // EmptyMass: XXX: StarSem says 6545 kg:
    constexpr static Mass   EmptyMass     = 6450.0_kg;

    // Mass of the RD-108A engine:
    constexpr static Mass   EngMass       = 1075.0_kg;

    // Masses of Fuel (Naftil) and Oxidiser:
    constexpr static Mass   FuelMass      = 26794.0_kg; // StarSem: 26300 (T1)
    constexpr static Mass   OxidMass      = 63709.0_kg; // StarSem: 63800
    constexpr static Mass   PerOxMass     = 2636.0_kg;
    constexpr static Mass   GasesMass     = 509.0_kg;
    constexpr static Mass   FullMass      = EmptyMass + FuelMass + OxidMass +
                                            PerOxMass + GasesMass;

    // UnSpendable Remnants of the Fuel and Oxidiser in Stage2 at the engine
    // cut-off time:
    constexpr static Mass   FuelRem       = 272.0_kg;
    constexpr static Mass   OxidRem       = 678.0_kg;
    constexpr static Mass   PerOxRem      = 263.0_kg;
    constexpr static Mass   GasesRem      = 84.0_kg;

    //-----------------------------------------------------------------------//
    // RD-108A (14D21) Engine Performance:                                   //
    //-----------------------------------------------------------------------//
    // Isp (SL/Vac, sec): EnergoMash says 257.7/320.6,  StarSem: 255.0/319.0;
    // but it is unclear whether these vals apply to the Main Engine    (w/o 
    // Vernier Chambers) or to the Whole Engine (incl Vernier Chambers, which
    // may have a slightly lower Isp). So assume higher vals for the Main En-
    // gine:
    constexpr static Time   IspMainSL     = 262.9_sec;
    constexpr static Time   IspMainVac    = 327.1_sec;

    // For Thrust (SL/Vac, tf), many different vals exist:
    // 70/87, 79.2/92.1, 80.8/94 (EnergoMash), 80.81/100.97;
    // but we must have IspVac/IspSL = ThrustVac/ThrustSL = 1.244 for the Main
    // Engine (see above);   then the 2nd and 3rd vals are clearly incorrect
    // (the ratio for them is 1.163,  so they are probably made of a mix  of
    // Main Engine and Vernier Chamber vals); the 1st and 4th vals yield reas-
    // onable ratios;  then the 1st one is apparently for Main Engine only,
    // and the last one is for the engine over-all;   use the former here:
    //
    constexpr static Force  ThrustMainSL  = 70000.0_kg * g0;
    constexpr static Force  ThrustMainVac =
      ThrustMainSL * double(IspMainVac / IspMainSL);  // 87.1 tf

    // For each Vernier Chamber, the data are taken from the old RD-107, so
    // they might not be very accurate. For example, the RD-107 Vernier Chamber
    // Fuel/Oxid Rates of 4.15/8.55 kg/sec would then lead to thrust of 3.2/4.0
    // tf, and the total Engine thrust (Main + 4 Vernier Chambers) of  82.8/103
    // tf, which is higher by 2 tf than the above figure of 80.81/100.97; so we
    // keep the Isp vals (their proportion is almost the same as for the  Main
    // Engine) and adjust the MassRate:
    // IMPORTANT:
    // Vernier Chambers of Stage2 are installed in positions at +-Y, +-Z  axes,
    // ie towards the corresp Stage1 blocks. They are rotatable to +-45 degress
    // in the planes parallel to the X axis, ie, the Vernier Chambers installed
    // at +-Y are rotatable in the XZ plane, and those installed at +-Z are ro-
    // tatable in the XY plane:
    //
    constexpr static Time   IspVernSL     = 251.9_sec;
    constexpr static Time   IspVernVac    = 313.1_sec;
    constexpr static Force  ThrustVernSL  = 2700.0_kg * g0;
    constexpr static Force  ThrustVernVac =
      ThrustVernSL * double(IspVernVac / IspVernSL);
      // 3.356 tf => Total=80.8/100.5 tf, very close to 80.81/100.97 tf above!

    // We can then calculate the MassRates for the Main Engine and for each
    // Vernier Chamber:
    constexpr static auto   MassRateMain  = ThrustMainSL  / (IspMainSL  * g0);
    static_assert(MassRateMain.ApproxEquals(ThrustMainVac / (IspMainVac * g0)));
    // ~277.89 kg/sec

    constexpr static auto   MassRateVern  = ThrustVernSL  / (IspVernSL  * g0);
    static_assert(MassRateVern.ApproxEquals(ThrustVernVac / (IspVernVac * g0)));
    // ~10.72  kg/sec, NOT 4.15+8.55=12.70 kg/sec as it would be with the orig-
    // inal RD-107 Vernier Chamber data...

    // Total MassRate for the whole Engine:
    constexpr static auto   MassRate      = MassRateMain + 4.0 * MassRateVern;
    // ~320.76 kg/sec

    // Separate Fuel and Oxid Rates are obtained using the Oxid/Fuel Ratio which
    // we derive from the over-all Fuel and Oxid spendable masses:
    constexpr static double OxidFuelRatio =
      double((OxidMass - OxidRem) / (FuelMass - FuelRem));
    // ~2.38, very close to 2.39 often quoted for RD-108A

    constexpr static auto   MassRateFuel  =
      MassRate                  / (1.0 + OxidFuelRatio);
    constexpr static auto   MassRateOxid  =
      MassRate * (OxidFuelRatio / (1.0 + OxidFuelRatio));

  private:
/*
    //-----------------------------------------------------------------------//
    // "MechElement"s: "Proto"s with Yet-UnKnown Masses:                     //
    //-----------------------------------------------------------------------//
    // XXX: Inter-Stage Grid (between Stages 2 and 3) is not modeled. Its mass
    // will be implicitly included into the EmptyMass of Stage2.
    // Then we need the X-Coord of Stage2 (somewhere) for reference. For that,
    // we take the Base (
    // Deflector Cone at the top of Stage2:
    constexpr static TrCone DeflectorCone =
      TrCone
      (
        ForeX0, 
        D,
        ForeH,
        Density(0.0)                    // No Propellant there
      );

  public:
    //-----------------------------------------------------------------------//
    // "MechElement"s with Real Masses:                                      //
    //-----------------------------------------------------------------------//
    // Aft Section (jettisonable): Cylinder:
    constexpr static TrCone AftSection =
      TrCone
      (
        OxidTankMidProto.GetRight()[0], // Common base with OxidTankMid (right)
        D,
        AftH,
        Density(0.0),                   // No Propellant in this Section
        AftMass                         // Mass is known!
      );

    // RD-0124 Engine, modeled as a PointMass:
    constexpr static PointMass Engine =
      PointMass
      (
        OxidTankMidProto.GetRight()[0] + EngCoMdX,
        0.0_m,
        0.0_m,
        EngMass                         // Mass is known!
      );

  private:
    // Scale Factor to determine the individual masses. NB: For simplicity, we
    // include the Mass of the Gases (which is very small) into the  TotalMass
    // of the "ConstElement"s:
    //
    constexpr static double ScaleFactor =
      ME::GetMassScale
      (
        // Components for which we provide the TotalMass (below):
        { &ForeSectionProto,
          &FuelTankUpProto, &FuelTankMidProto, &FuelTankLowProto,
          &EquipBayProto,
          &OxidTankUpProto, &OxidTankMidProto, &OxidTankLowProto
        },
        // The TotalMass (incl Gases):
        EmptyMass + GasesMass - AftMass - EngMass
      );

  public:
    // Real-Mass Fore Section:
    constexpr static TrCone    ForeSection =
      ME::ProRateMass(ForeSectionProto, ScaleFactor);

    // Real-Mass Fuel Tank Components:
    constexpr static SpherSegm FuelTankUp  =
      ME::ProRateMass(FuelTankUpProto,  ScaleFactor);

    constexpr static TrCone    FuelTankMid =
      ME::ProRateMass(FuelTankMidProto, ScaleFactor);

    constexpr static SpherSegm FuelTankLow =
      ME::ProRateMass(FuelTankLowProto, ScaleFactor);

    // Real-Mass Equipment Bay (XXX: though the mass  of  the Control System
    // itself cannot be exactly accounted for -- it is pro-rated in the same
    // way as the masses of all other components):
    constexpr static TrCone    EquipBay    =
      ME::ProRateMass(EquipBayProto,    ScaleFactor);

    // Real-Mass Oxid Tank Components:
    constexpr static SpherSegm OxidTankUp  =
      ME::ProRateMass(OxidTankUpProto,  ScaleFactor);

    constexpr static TrCone    OxidTankMid =
      ME::ProRateMass(OxidTankMidProto, ScaleFactor);

    constexpr static SpherSegm OxidTankLow =
      ME::ProRateMass(OxidTankLowProto, ScaleFactor);

    // "MechElement"s for the Empty Stage with Gases ("EG"): Before and After
    // Jettisoning of the AftSection:
    constexpr static ME EGAfterME  =
      ForeSection + FuelTankUp  + FuelTankMid  + FuelTankLow +
      EquipBay    + OxidTankUp  + OxidTankMid  + OxidTankLow + Engine;

    constexpr static ME EGBeforeME = EGAfterME + AftSection;
 
    //-----------------------------------------------------------------------//
    // Geometry Checks:                                                      //
    //-----------------------------------------------------------------------//
    // The bottom of the FuelTank must not touch the top of the Oxidiser Tank;
    // there must be a positive gap between them:
    //
    constexpr static Len FuelOxidTanksGap =
      OxidTankUp.GetLeft()[0] - FuelTankLow.GetRight()[0];
    static_assert(IsPos (FuelOxidTanksGap));

    // On the other hand, the top of the FuelTank appears over the over-all
    // cylindrical shell (but inside both kinds of InterStages):
    constexpr static Len FuelTankTop  = FuelTankUp.GetLeft()[0];
    static_assert(IsNeg (FuelTankTop) &&
                  FuelTankTop > - SH::InterStageLargeH  &&
                  FuelTankTop > - SH::InterStageSmallH);

    // Also, the bottom of the OxidiserTank must be well inside the over-all
    // cylindrical shell:
    static_assert(OxidTankLow.GetRight()[0] < H);

  private:
    //-----------------------------------------------------------------------//
    // "MechElement"s for Propellant Loads:                                  //
    //-----------------------------------------------------------------------//
    // Propallant Mass Capacities (MC) of Fuel and Oxid Tank Sections:
    constexpr static Mass FuelTankUpMC     = FuelTankUp .GetPropMassCap();
    constexpr static Mass FuelTankMidMC    = FuelTankMid.GetPropMassCap();
    constexpr static Mass FuelTankLowMC    = FuelTankLow.GetPropMassCap();
    constexpr static Mass FuelTankLowMidMC = FuelTankLowMC + FuelTankMidMC;

    constexpr static Mass OxidTankUpMC     = OxidTankUp .GetPropMassCap();
    constexpr static Mass OxidTankMidMC    = OxidTankMid.GetPropMassCap();
    constexpr static Mass OxidTankLowMC    = OxidTankLow.GetPropMassCap();
    constexpr static Mass OxidTankLowMidMC = OxidTankLowMC + OxidTankMidMC;

  public:
    // Volumes of the Fuel and Oxid Tanks:
    constexpr static Vol  FuelTankVol  =
      FuelTankUp .GetEnclVol() + FuelTankMid.GetEnclVol() +
      FuelTankLow.GetEnclVol();

    constexpr static Vol  OxidTankVol  =
      OxidTankUp .GetEnclVol() + OxidTankMid.GetEnclVol() +
      OxidTankLow.GetEnclVol();

    // Maximum Theoretical Fuel and Oxid Capacities of the resp Tanks:
    constexpr static Mass FuelTankMC   =
      FuelTankUpMC + FuelTankMidMC + FuelTankLowMC;
    constexpr static Mass OxidTankMC   =
      OxidTankUpMC + OxidTankMidMC + OxidTankLowMC;

  private:
    // "ConstElement" objs for Maximum Theoretical Propellant Loads in Tank
    // Sections (for optimisation of "GetDynParams"):
    constexpr static ME   FuelUpME     = FuelTankUp .GetPropBulkME(FuelTankUpMC);
    constexpr static ME   FuelMidME    = FuelTankMid.GetPropBulkME(FuelTankMidMC);
    constexpr static ME   FuelLowME    = FuelTankLow.GetPropBulkME(FuelTankLowMC);
    constexpr static ME   FuelLowMidME = FuelLowME + FuelMidME;

    constexpr static ME   OxidUpME      = OxidTankUp .GetPropBulkME(OxidTankUpMC);
    constexpr static ME   OxidMidME     = OxidTankMid.GetPropBulkME(OxidTankMidMC);
    constexpr static ME   OxidLowME     = OxidTankLow.GetPropBulkME(OxidTankLowMC);
    constexpr static ME   OxidLowMidME  = OxidLowME + OxidMidME;

  public:
    // Fuel and Oxid Load Ratios (ActualMass / TheorMassCapacity):
    constexpr static double FuelLoadRatio = double(FuelMass / FuelTankMC);
    static_assert(FuelLoadRatio < 1.0);

    constexpr static double OxidLoadRatio = double(OxidMass / OxidTankMC);
    static_assert(OxidLoadRatio < 1.0);

    //-----------------------------------------------------------------------//
    // TimeLine Consts:                                                      //
    //-----------------------------------------------------------------------//
    // Once again, "FT" stands for "FullThrust":
    //
    constexpr static Time IgnTime     = SC::Stage2IgnTime;
    constexpr static Time FTStartTime = SC::Stage2FullThrustTime;
    constexpr static Time AftJetTime  = SC::Stage2AftJetTime;
    constexpr static Time FTEndTime   = FTStartTime + SC::Stage2FTBurnDur;
    constexpr static Time CutOffTime  = FTEndTime   + SC::Stage2ThrBurnDur;

    // Ordering of Events:
    static_assert
      (IsPos(IgnTime)           && IgnTime    < FTStartTime &&
       FTStartTime < AftJetTime && AftJetTime < FTEndTime   &&
       FTEndTime   < CutOffTime);

    //-----------------------------------------------------------------------//
    // Max Theoretical Burn Duration at Full Thrust (FT):                    //
    //-----------------------------------------------------------------------//
    // NB: The following takes into account 2 periods of Throttled run as well,
    // which reduce the FullThrust Burn Duration:
    //
    constexpr static Time   MaxFTBurnDur =
      std::min((FuelMass - FuelRem) / FuelRateFT,
               (OxidMass - OxidRem) / OxidRateFT)
      - 2.0 * SC::Stage2ThrBurnDur * ThrottlRate;

    // The actual FullThrust Burn Duration must not exceed the above theoreti-
    // cal maximum:
    constexpr static Time   FTBurnDur = SC::Stage2FTBurnDur;
    static_assert(FTBurnDur < MaxFTBurnDur);

  private:
    //-----------------------------------------------------------------------//
    // Masses at Event Times (for "GetDynParams" Optimisation):              //
    //-----------------------------------------------------------------------//
    constexpr static auto ThrMassRate = ThrottlRate * MassRateFT;
    constexpr static auto ThrFuelRate = ThrottlRate * FuelRateFT;
    constexpr static auto ThrOxidRate = ThrottlRate * OxidRateFT;

    // Masses at FullThrust Start: Full, Fuel and Oxid:
    constexpr static Mass FTStartFullMass =
      FullMass        - ThrMassRate * (SC::Stage2FullThrustTime  - IgnTime);

    constexpr static Mass FTStartFuelMass =
      FuelMass        - ThrFuelRate * (SC::Stage2FullThrustTime  - IgnTime);

    constexpr static Mass FTStartOxidMass =
      OxidMass        - ThrOxidRate * (SC::Stage2FullThrustTime  - IgnTime);

    // Full Mass just after jettisoning of the Aft Section:
    constexpr static Mass AftJetFullMass  =
      FTStartFullMass - MassRateFT * (AftJetTime - FTStartTime) - AftMass;

    // Masses at the end of Full Thrust: Full, Fuel and Oxid:
    constexpr static Mass FTEndFullMass   =
      AftJetFullMass  - MassRateFT * (FTEndTime  - AftJetTime);

    constexpr static Mass FTEndFuelMass   =
      FTStartFuelMass - FuelRateFT * (FTEndTime  - FTStartTime);

    constexpr static Mass FTEndOxidMass   =
      FTStartOxidMass - OxidRateFT * (FTEndTime  - FTStartTime);

  public:
    // Masses at the Engine Cut-Off:
    constexpr static Mass CutOffFullMass  =
      FTEndFullMass   - ThrMassRate  * SC::Stage2ThrBurnDur;

    constexpr static Mass CutOffFuelMass  =
      FTEndFuelMass   - ThrFuelRate  * SC::Stage2ThrBurnDur;

    constexpr static Mass CutOffOxidMass  =
      FTEndOxidMass   - ThrOxidRate  * SC::Stage2ThrBurnDur;

  private:
    // Checks: The Empty incl Gases ("RG") Mass must be the same in all cases,
    // up to the mass of the AftSection (jettisonable):
    //
    constexpr static Mass EGMassBefore  = EmptyMass    + GasesMass;
    constexpr static Mass EGMassAfter   = EGMassBefore - AftMass;
    static_assert  (IsPos(EGMassAfter));

    constexpr static Mass FTStartEGMass =
      FTStartFullMass - FTStartFuelMass - FTStartOxidMass;

    constexpr static Mass FTEndEGMass   =
      FTEndFullMass   - FTEndFuelMass   - FTEndOxidMass;

    constexpr static Mass CutOffEGMass  =
      CutOffFullMass  - CutOffFuelMass  - CutOffOxidMass;

    static_assert(FTStartEGMass.ApproxEquals(EGMassBefore) &&
                  FTEndEGMass  .ApproxEquals(EGMassAfter)  &&
                  CutOffEGMass .ApproxEquals(EGMassAfter));

    // Also, "EGMass{Before,After}" must be equal to the corresp vals in the
    // "EG{Before,After}ME":
    static_assert(EGMassBefore.ApproxEquals(EGBeforeME.GetMass()) &&
                  EGMassAfter .ApproxEquals(EGAfterME .GetMass()));

    // Also, Fuel and Oxid Masses at CutOff must be greater than UnSpendable
    // Remnants:
    static_assert(CutOffFuelMass > FuelRem &&
                  CutOffOxidMass > OxidRem);

  public:
    //=======================================================================//
    // Dynamic Params as a function of Flight Time:                          //
    //=======================================================================//
    // NB: This method is NOT "constexpr": it is intended to be called at Run-
    // Time (eg multiple times during Trajectory Integration):
    //
    static StageDynParams GetDynParams(Time a_t);
*/
  };
}
