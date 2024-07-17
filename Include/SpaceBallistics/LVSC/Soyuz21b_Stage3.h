// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/LVSC/Soyuz21b_Stage3.h":                 //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/TrConeSpherSegm.hpp"
#include "SpaceBallistics/LVSC/LVSC.h"
#include "SpaceBallistics/LVSC/Soyuz21b_Consts.h"
#include "SpaceBallistics/LVSC/Soyuz21b_Head.h"
#include "SpaceBallistics/LVSC/Propellants.h"
#include "SpaceBallistics/LVSC/StageDynParams.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC  = Soyuz21b_Consts;
  namespace SH  = Soyuz21b_Head;

  // All "MechElements" are instantiated with "LVSC::Soyuz21b":
  using     ME  = MechElement<LVSC::Soyuz21b>;
  using     PM  = PointMass  <LVSC::Soyuz21b>;
  using     TC  = TrCone     <LVSC::Soyuz21b>;
  using     SpS = SpherSegm  <LVSC::Soyuz21b>;

  //=========================================================================//
  // "Soyuz21b_Stage3" Class:                                                //
  //=========================================================================//
  // The following local co-ord system is used for convenience (so that most
  // X-coords are positive and the origin  is not  affected by separation of
  // stages):
  // (*) The OX axis is the main axis of the rocket. The OX positive direction
  //     is towards the HEAD of the LV. The origin O is the LOWER base of the
  //     InterStage (junction plane with Stage3). That is, X >= 0 for
  //     for InterStage, the optional Stage4, Payload Adapter/Dispenser, Pay-
  //     load itself and the Fairing, and X <= 0 for Stages 3, 2, 1.
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
    constexpr static Len  H               = SC::Stage3Len;

    // Stage3 Fore Section:
    constexpr static Len  ForeX0          = SC::X0;    // THE ORIGIN IS HERE!!!
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
    // resp., but those figures seem to be GROSSLY INCORRECT (too high for the
    // actual tank volumes!):
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

    // IMPORTANT:
    // Each of the 4 Chambers of RD-0124 is GIMBALED in the Tangential  Plane.
    // The Chambers are installed in the same planes as the Blocks B,V,G,D of
    // Stage1, that is, in planes XY and XZ, at some distance (r > 0) from the
    // Stage 3 main axis (the X axis):
    //    Chamber0: (Y= r, Z= 0);
    //    Chamber1: (Y= 0, Z= r);
    //    Chamber2: (Y=-r, Z= 0);
    //    Chamber3: (Y= 0, Z=-r).
    // The max gimbaling angle of RD-0124 Chambers is not known exactlt, prob-
    // ably in the range +-(3..8) degs; we can safely assume +-5 degs.
    // BY OUR CONVENTION, GimbalAngle > 0 corresponds to a COUNTER-CLOCK-WISE
    // rotation of the corresp Chamber (more precisely, NOZZLE) in the YZ plane
    // around the X axis.
    // NORMALLY, one can expect that Chambers (0 and 2), (1 and 3) are moved
    // symmetrically in one direction, so
    //    GimbalAngle[0] = - GimbalAngle[2],
    //    GimbalAngle[1] = - GimbalAngle[3],
    // but this is not strictly enforced:
    //
    constexpr static Angle_deg GimbalAmpl = Angle_deg(5.0);
    using            ThrustVecCtl         = std::array<Angle_deg, 4>;

  private:
    //-----------------------------------------------------------------------//
    // "MechElement"s: "Proto"s with Yet-UnKnown Masses:                     //
    //-----------------------------------------------------------------------//
    // Fore Section:
    constexpr static TC ForeSectionProto  =
      TC
      (
        ForeX0,                         // Up @ X=0
        D,
        ForeH,
        Density(0.0)                    // No Propellant there
      );

    // Fuel Tank "Proto" components (with yet-unknown masses):
    constexpr static SpS FuelTankUpProto  =
      SpS
      (
        true,                           // Facing Up
        ForeSectionProto.GetLow()[0],   // Base @ Low (Bottom) of ForeSection
        D,
        Propellants::RG1Dens            // Naftil
      );
    constexpr static TC  FuelTankMidProto =
      TC
      (
        ForeSectionProto.GetLow()[0],   // Up   @ Low (Bottom) of ForeSection
        D,
        FuelTankMidH,
        Propellants::RG1Dens            // Naftil
      );
    constexpr static SpS FuelTankLowProto =
      SpS
      (
        false,                          // Facing Low
        FuelTankMidProto.GetLow()[0],   // Base @ Low (Bottom) of FuelTankMid
        D,
        FuelTankLowH,
        Propellants::RG1Dens            // Naftil
      );

    // Equipment Bay (containing the Control System):
    // "Proto" with yet-unknown mass:
    constexpr static TC  EquipBayProto =
      TC
      (
        FuelTankMidProto.GetLow()[0],   // Up @ Low (Bottom) of FuelTankMid
        D,
        EquipBayH,
        Density(0.0)                    // No Propellant there
      );

    // Oxidiser Tank "Proto" components (with yet-unknown masses):
    constexpr static SpS OxidTankUpProto =
      SpS
      (
        true,                           // Facing Up
        EquipBayProto.GetLow()[0],      // Base @ Low (Bottom) of EquipBay
        D,
        Propellants::LOxDens
      );
    constexpr static TC  OxidTankMidProto =
      TC
      (
        EquipBayProto.GetLow()[0],      // Up @   Low (Bottom) of EquipBay
        D,
        OxidTankMidH,
        Propellants::LOxDens
      );
    constexpr static SpS OxidTankLowProto =
      SpS
      (
        false,                          // Facing Low
        OxidTankMidProto.GetLow()[0],   // Base @ Low (Bottom) of OxidTankMid
        D,
        Propellants::LOxDens
      );

  public:
    //-----------------------------------------------------------------------//
    // "MechElement"s with Real Masses:                                      //
    //-----------------------------------------------------------------------//
    // Aft Section (jettisonable): Cylinder:
    constexpr static TC  AftSection =
      TC
      (
        OxidTankMidProto.GetLow()[0],   // Base @ Low (Bottom) of OxidTankMid
        D,
        AftH,
        Density(0.0),                   // No Propellant in this Section
        AftMass                         // Mass is known!
      );

    // RD-0124 Engine, modeled as a PointMass:
    constexpr static PM  Engine =
      PM
      (
        OxidTankMidProto.GetLow()[0] - EngCoMdX,
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
    constexpr static TC  ForeSection =
      ME::ProRateMass(ForeSectionProto, ScaleFactor);

    // Real-Mass Fuel Tank Components:
    constexpr static SpS FuelTankUp  =
      ME::ProRateMass(FuelTankUpProto,  ScaleFactor);

    constexpr static TC  FuelTankMid =
      ME::ProRateMass(FuelTankMidProto, ScaleFactor);

    constexpr static SpS FuelTankLow =
      ME::ProRateMass(FuelTankLowProto, ScaleFactor);

    // Real-Mass Equipment Bay (XXX: though the mass  of  the Control System
    // itself cannot be exactly accounted for -- it is pro-rated in the same
    // way as the masses of all other components):
    constexpr static TC  EquipBay    =
      ME::ProRateMass(EquipBayProto,    ScaleFactor);

    // Real-Mass Oxid Tank Components:
    constexpr static SpS OxidTankUp  =
      ME::ProRateMass(OxidTankUpProto,  ScaleFactor);

    constexpr static TC  OxidTankMid =
      ME::ProRateMass(OxidTankMidProto, ScaleFactor);

    constexpr static SpS OxidTankLow =
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
      FuelTankLow.GetLow()[0] - OxidTankUp.GetUp()[0];
    static_assert(IsPos (FuelOxidTanksGap));

    // On the other hand, the top of the FuelTank appears over the over-all
    // cylindrical shell (but inside both kinds of InterStages):
    constexpr static Len FuelTankTop  = FuelTankUp.GetUp()[0];
    static_assert(IsPos (FuelTankTop) &&
                  FuelTankTop < ForeX0 + SH::InterStageLargeH  &&
                  FuelTankTop < ForeX0 + SH::InterStageSmallH);

    // Also, the bottom of the OxidiserTank must be well inside the over-all
    // cylindrical shell:
    static_assert(OxidTankLow.GetLow()[0] > ForeX0 - H);

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
    constexpr static ME FuelUpME     = FuelTankUp .GetPropBulkME(FuelTankUpMC) ;
    constexpr static ME FuelMidME    = FuelTankMid.GetPropBulkME(FuelTankMidMC);
    constexpr static ME FuelLowME    = FuelTankLow.GetPropBulkME(FuelTankLowMC);
    constexpr static ME FuelLowMidME = FuelLowME + FuelMidME;

    constexpr static ME OxidUpME     = OxidTankUp .GetPropBulkME(OxidTankUpMC) ;
    constexpr static ME OxidMidME    = OxidTankMid.GetPropBulkME(OxidTankMidMC);
    constexpr static ME OxidLowME    = OxidTankLow.GetPropBulkME(OxidTankLowMC);
    constexpr static ME OxidLowMidME = OxidLowME + OxidMidME;

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
    static StageDynParams<LVSC::Soyuz21b>
    GetDynParams(Time a_t, ThrustVecCtl const& a_thrust_ctl);
  };
}
// End namespace SpaceBallistics
