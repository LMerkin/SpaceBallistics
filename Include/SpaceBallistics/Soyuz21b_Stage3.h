// vim:ts=2:et
//===========================================================================//
//                             "Soyuz21b_Stage3.h":                          //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ConstrElement.hpp"
#include "SpaceBallistics/Soyuz21b_Consts.h"
#include "SpaceBallistics/Soyuz21b_Head.h"
#include "SpaceBallistics/Propellants.h"
#include "SpaceBallistics/StageDynParams.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC = Soyuz21b_Consts;
  namespace SH = Soyuz21b_Head;
  using     CE = ConstrElement;

  //=========================================================================//
  // "Soyuz21b_Stage3" Class:                                                //
  //=========================================================================//
  // The following local co-ord system is used for convenience (so that most
  // X-coords are positive and the origin  is not  affected by separation of
  // stages):
  // (*) The OX axis is the main axis of the rocket. The OX positive direction
  //     is towards the TAIL (REAL of the LV). The origin O is the LOWER base
  //     of the InterStage (junction plane with Stage3).  That is, X <= 0 for
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
  class Soyuz21b_Stage3
  {
  private:
    //=======================================================================//
    // No objects can be created of this class:                              //
    //=======================================================================//
    Soyuz21b_Stage3() = delete;

  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Geometry:                                                             //
    //-----------------------------------------------------------------------//
    // Stage3 Over-All:
    constexpr static Len  D               = SH::InterStageD1;
    constexpr static Len  H               = 6.745_m;

    // Stage3 Fore Section:
    constexpr static Len  ForeX0          = 0.0_m;  // THE ORIGIN IS HERE!!!
    constexpr static Len  ForeH           = 0.718_m;

    // Stage3 Fuel Tank:
    constexpr static Len  FuelTankMidH    = 0.373_m;
    constexpr static Len  FuelTankLowH    = 0.552_m;

    // Stage3 Equipment Bay:
    constexpr static Len  EquipBayH       = 2.071_m;

    // Stage3 Oxidiser Tank:
    constexpr static Len  OxidTankMidH    = 1.05_m;

    // Stage3 Aft Section:
    constexpr static Len  AftH            =
      H - (ForeH + FuelTankMidH + EquipBayH + OxidTankMidH);
    static_assert(IsPos(AftH));

    // XXX: The approximate dX of the CoM of the Engine from AftTop:
    constexpr static Len  EngCoMdX        = 1.0_m;

    //-----------------------------------------------------------------------//
    // Masses:                                                               //
    //-----------------------------------------------------------------------//
    // EmptyMass (incl the jettisonable Aft Section and the RD-0124 engine):
    // XXX: Some srcs give a much larger EmptyMass: 2597..2710 kg, or a lower
    // one: 2355 kg:
    constexpr static Mass   EmptyMass     = 2490.0_kg;

    // Mass of the jettisonable Aft Section:
    constexpr static Mass   AftMass       = 441.0_kg;

    // Mass of the RD-0124 engine:
    // Data in various srcs vary significantly: 450 kg, 460 kg, 520 kg, 572 kg;
    // we assume that the latter figure includes the support structure:
    constexpr static Mass   EngMass       = 572.0_kg;

    // Masses of Fuel and Oxidiser. ArianSpace/StarSem says 7600 and 17800 kg,
    // resp., but those figures seem to be incorrect  (too much for the actual
    // tank volumes!):
    constexpr static Mass   FuelMass      = 6650.0_kg;
    constexpr static Mass   OxidMass      = 16554.0_kg;
    constexpr static Mass   GasesMass     = 30.0_kg;  // He, air, ...
    constexpr static Mass   FullMass      = EmptyMass + FuelMass + OxidMass +
                                            GasesMass;

    // UnSpendable Remnants of the Fuel and Oxidiser in Stage3 at the engine
    // cut-off time:
    constexpr static Mass   FuelRem       = 104.0_kg; // StarSem: 98 kg
    constexpr static Mass   OxidRem       = 167.0_kg; // StarSem: 188..207 kg

    //-----------------------------------------------------------------------//
    // RD-0124 (14D23) Engine Performance:                                   //
    //-----------------------------------------------------------------------//
    // (In Vacuum; SeaLevel is meaningless here). StarSem: Thrust = 297.9e3 N
    constexpr static Time   IspVac        = 359.0_sec;
    constexpr static Force  ThrustVac     = Force(294.3e3);

    // The actual MassRates:
    // XXX: These figures may not be highly accurate, as it is unclear whether
    // they refer to Naftil or Kerosene. The OxidRate is quoted somewhere as
    // 56.7 kg/sec, which is clearly too low.
    // In different sources, Oxidiser/Fuel Ratio is 2.5..2.6, here 2.50:
    // These figure must be consistent with the BurnDur and the StaticThrust
    // below:
    constexpr static auto   FuelRate      = 23.8_kg / 1.0_sec;
    constexpr static auto   OxidRate      = 59.6_kg / 1.0_sec;
    constexpr static auto   MassRate      = FuelRate + OxidRate;

    // The MassRate is connected to Specific Impulse and Thrust, but we must
    // take into account that Thrust is a sum of Reactive Force (proportional
    // to MassRate) and the residual inline pressure force at the nozzle exh-
    // aust:
    constexpr static Force  StaticThrust  =
      ThrustVac - MassRate * IspVac * g0;
    static_assert(IsPos(StaticThrust));

  private:
    //-----------------------------------------------------------------------//
    // "ConstrElement"s: "Proto"s with Yet-UnKnown Masses:                   //
    //-----------------------------------------------------------------------//
    // Fore Section:
    constexpr static TrCone ForeSectionProto =
      TrCone
      (
        ForeX0,
        D,
        ForeH,
        Density(0.0)                    // No Propellant there
      );

    // Fuel Tank "Proto" components (with yet-unknown masses):
    constexpr static SpherSegm FuelTankUpProto  =
      SpherSegm
      (
        false,                          // Facing Up = Left = !Right
        ForeSectionProto.GetRight()[0], // Bottom of ForeSection
        D,
        RG1Dens                         // Naftil
      );
    constexpr static TrCone    FuelTankMidProto =
      TrCone
      (
        ForeSectionProto.GetRight()[0], // Common base with FuelTankUp
        D,
        FuelTankMidH,
        RG1Dens                         // Naftil
      );
    constexpr static SpherSegm FuelTankLowProto =
      SpherSegm
      (
        true,                           // Facing Down = Right
        FuelTankMidProto.GetRight()[0], // Common base with FuelTankMid (right)
        D,
        FuelTankLowH,
        RG1Dens                         // Naftil
      );

    // Equipment Bay (containing the Control System):
    // "Proto" with yet-unknown mass:
    constexpr static TrCone    EquipBayProto =
      TrCone
      (
        FuelTankMidProto.GetRight()[0], // Common base with FuelTankMid (right)
        D,
        EquipBayH,
        Density(0.0)                    // No Propellant there
      );

    // Oxidiser Tank "Proto" components (with yet-unknown masses):
    constexpr static SpherSegm OxidTankUpProto =
      SpherSegm
      (
        false,                          // Facing Up = Left = !Right
        EquipBayProto.GetRight()[0],    // Common base with EquipBay (right)
        D,
        LOxDens
      );
    constexpr static TrCone    OxidTankMidProto =
      TrCone
      (
        EquipBayProto.GetRight()[0],    // Common base with EquipBay (right)
        D,
        OxidTankMidH,
        LOxDens
      );
    constexpr static SpherSegm OxidTankLowProto =
      SpherSegm
      (
        true,                           // Facing Down = Right
        OxidTankMidProto.GetRight()[0], // Common base with OxidTankMid (right)
        D,
        LOxDens
      );

  public:
    //-----------------------------------------------------------------------//
    // "ConstrElement"s with Real Masses:                                    //
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
      CE::GetMassScale
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
      CE::ProRateMass(ForeSectionProto, ScaleFactor);

    // Real-Mass Fuel Tank Components:
    constexpr static SpherSegm FuelTankUp  =
      CE::ProRateMass(FuelTankUpProto,  ScaleFactor);

    constexpr static TrCone    FuelTankMid =
      CE::ProRateMass(FuelTankMidProto, ScaleFactor);

    constexpr static SpherSegm FuelTankLow =
      CE::ProRateMass(FuelTankLowProto, ScaleFactor);

    // Real-Mass Equipment Bay (XXX: though the mass  of  the Control System
    // itself cannot be exactly accounted for -- it is pro-rated in the same
    // way as the masses of all other components):
    constexpr static TrCone    EquipBay    =
      CE::ProRateMass(EquipBayProto,    ScaleFactor);

    // Real-Mass Oxid Tank Components:
    constexpr static SpherSegm OxidTankUp  =
      CE::ProRateMass(OxidTankUpProto,  ScaleFactor);

    constexpr static TrCone    OxidTankMid =
      CE::ProRateMass(OxidTankMidProto, ScaleFactor);

    constexpr static SpherSegm OxidTankLow =
      CE::ProRateMass(OxidTankLowProto, ScaleFactor);

    // "ConstrElement"s for the Empty Stage with Gases ("EG"): Before and After
    // Jettisoning of the AftSection:
    constexpr static CE EGAfterCE  =
      ForeSection + FuelTankUp  + FuelTankMid  + FuelTankLow +
      EquipBay    + OxidTankUp  + OxidTankMid  + OxidTankLow + Engine;

    constexpr static CE EGBeforeCE = EGAfterCE + AftSection;

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
    // "ConstrElement"s for Propellant Loads:                                //
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
    constexpr static CE   FuelUpCE     = FuelTankUp .GetPropCE(FuelTankUpMC) ;
    constexpr static CE   FuelMidCE    = FuelTankMid.GetPropCE(FuelTankMidMC);
    constexpr static CE   FuelLowCE    = FuelTankLow.GetPropCE(FuelTankLowMC);
    constexpr static CE   FuelLowMidCE = FuelLowCE + FuelMidCE;

    constexpr static CE   OxidUpCE     = OxidTankUp .GetPropCE(OxidTankUpMC) ;
    constexpr static CE   OxidMidCE    = OxidTankMid.GetPropCE(OxidTankMidMC);
    constexpr static CE   OxidLowCE    = OxidTankLow.GetPropCE(OxidTankLowMC);
    constexpr static CE   OxidLowMidCE = OxidLowCE + OxidMidCE;

  public:
    // Fuel and Oxid Load Ratios (ActualMass / TheorMassCapacity):
    constexpr static double FuelLoadRatio = double(FuelMass / FuelTankMC);
    static_assert(FuelLoadRatio < 1.0);

    constexpr static double OxidLoadRatio = double(OxidMass / OxidTankMC);
    static_assert(OxidLoadRatio < 1.0);

    //-----------------------------------------------------------------------//
    // TimeLine Consts:                                                      //
    //-----------------------------------------------------------------------//
    constexpr static Time IgnTime     = SC::Stage3IgnTime;
    constexpr static Time AftJetTime  = SC::Stage3AftJetTime;
    constexpr static Time CutOffTime  = SC::Stage3CutOffTime;

    // Ordering of Events:
    static_assert
      (IsPos(IgnTime) && IgnTime < AftJetTime && AftJetTime < CutOffTime);

    //-----------------------------------------------------------------------//
    // Max Theoretical Burn Duration:                                        //
    //-----------------------------------------------------------------------//
    constexpr static Time MaxBurnDur =
      std::min((FuelMass - FuelRem) / FuelRate,
               (OxidMass - OxidRem) / OxidRate);

    // The actual Burn Duration must not exceed the above theoretical maximum:
    constexpr static Time BurnDur = SC::Stage3BurnDur;
    static_assert(BurnDur < MaxBurnDur);

  public:
    //-----------------------------------------------------------------------//
    // Masses at Event Times (for "GetDynParams" Optimisation):              //
    //-----------------------------------------------------------------------//
    // Full Mass just after jettisoning of the Aft Section:
    constexpr static Mass AftJetFullMass  =
      FullMass - MassRate * (AftJetTime - IgnTime) - AftMass;

    // Masses at the Engine Cut-Off:
    constexpr static Mass CutOffFullMass  =
      AftJetFullMass - MassRate  * (CutOffTime - AftJetTime);

    constexpr static Mass CutOffFuelMass  =
      FuelMass       - FuelRate  * SC::Stage3BurnDur;

    constexpr static Mass CutOffOxidMass  =
      OxidMass       - OxidRate  * SC::Stage3BurnDur;

  private:
    // Checks: The Empty incl Gases ("EG") Mass must be the same in all cases,
    // up to the mass of the AftSection (jettisonable):
    //
    constexpr static Mass EGMassBefore  = EmptyMass    + GasesMass;
    constexpr static Mass EGMassAfter   = EGMassBefore - AftMass;
    static_assert  (IsPos(EGMassAfter));

    constexpr static Mass IgnEGMass     = FullMass - FuelMass - OxidMass;

    constexpr static Mass CutOffEGMass  =
      CutOffFullMass - CutOffFuelMass   - CutOffOxidMass;

    static_assert(IgnEGMass   .ApproxEquals(EGMassBefore) &&
                  CutOffEGMass.ApproxEquals(EGMassAfter));

    // Also, "EGMass{Before,After}" must be equal to the corresp vals in the
    // "EG{Before,After}CE":
    static_assert(EGMassBefore.ApproxEquals(EGBeforeCE.GetMass()) &&
                  EGMassAfter .ApproxEquals(EGAfterCE .GetMass()));

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
  };
}
