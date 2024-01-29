// vim:ts=2:et
//===========================================================================//
//                            "Soyuz21b_Stage3.h":                           //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ConstrElement.hpp"
#include "SpaceBallistics/Soyuz21b_Consts.h"
#include "SpaceBallistics/Soyuz21b_Stage3_Consts.h"
#include "SpaceBallistics/StageDynParams.h"

namespace SpaceBallistics
{
  // Consts pertinent to the LV as a whole:
  namespace SC  = Soyuz21b_Consts;

  // Consts pertinent to Stage3:
  namespace S3C = Soyuz21b_Stage3_Consts;

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
  //
  class Soyuz21b_Stage3
  {
  private:

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Construction Elements:                                                //
    //-----------------------------------------------------------------------//
    // Fore Section:
    TrCone          m_foreSection;

    // Fuel Tank:
    SpherSegm       m_fuelTankUp;
    TrCone          m_fuelTankMid;
    SpherSegm       m_fuelTankLow;

    // Equipment Bay (containing the Control System):
    TrCone          m_equipBay;

    // Oxidiser Tank:
    SpherSegm       m_oxidTankUp;
    TrCone          m_oxidTankMid;
    SpherSegm       m_oxidTankLow;

  public:
    // Aft Section (jettisonable):
    TrCone    const m_aftSection;

    // RD-0124 Engine, modeled as a PointMass:
    PointMass const m_engine;

    // For Information Only:
    // Max Theoretical Volumes and Masses of Fuel and Oxidiser:
    Mass      const m_maxFuelMass;
    Mass      const m_maxOxidMass;
    // Max Theoretical Burn Time at Full Thrust:
    Time      const m_maxBurnTimeFT;

    //-----------------------------------------------------------------------//
    // Flight Profile Params:                                                //
    //-----------------------------------------------------------------------//
    // TimeLine of Stage3 Operation:
    constexpr static Time IgnTime      = SC::Stage3IgnTime;
    constexpr static Time FTStartTime  = SC::Stage3FullThrustTime;
    constexpr static Time AftJetTime   = SC::Stage3AftJetTime;
    Time      const  m_ftEndTime       = FTStartTime + SC::Stage3FTBurnDur;
    Time      const  m_cutOffTime      = m_ftEndTime + SC::Stage3ThrBurnDur;

    // Masses at Time Nodes:
    constexpr static Mass FullMass     = S3C::FullMass;
    constexpr static auto ThrMassRate  = S3C::Throttling  * S3C::MassRateFT;
    constexpr static auto ThrFuelRate  = S3C::Throttling  * S3C::FuelRateFT;
    constexpr static auto ThrOxidRate  = S3C::Throttling  * S3C::OxidRateFT;

    // Masses at Full Thrust Start instant: Full, Fuel and Oxid:
    constexpr static Mass FTStartMass  =
      FullMass - ThrMassRate * (SC::Stage3FullThrustTime - IgnTime);

    constexpr static Mass FTStartFuel  =
      S3C::FuelMass
               - ThrFuelRate * (SC::Stage3FullThrustTime - IgnTime);

    constexpr static Mass FTStartOxid  =
      S3C::OxidMass
               - ThrOxidRate * (SC::Stage3FullThrustTime - IgnTime);

    // Mass just after jettisoning of the Aft Section:
    constexpr static Mass AftJetMass   =
    FTStartMass -
      S3C::MassRateFT * (AftJetTime - FTStartTime) -
      S3C::AftMass;

    // Masses at the end of Full Thrust:
    Mass      const  m_ftEndMass       =
      AftJetMass  - S3C::MassRateFT * (m_ftEndTime - AftJetTime);

    Mass      const  m_ftEndFuel       =
      FTStartFuel - S3C::FuelRateFT * (m_ftEndTime - FTStartTime);

    Mass      const  m_ftEndOxid       =
      FTStartOxid - S3C::OxidRateFT * (m_ftEndTime - FTStartTime);

    // Masses at the Engine Cut-Off:
    Mass      const  m_cutOffMass      =
      m_ftEndMass - ThrMassRate    * SC::Stage3ThrBurnDur;

    Mass      const m_cutOffFuel       =
      m_ftEndFuel - ThrFuelRate    * SC::Stage3ThrBurnDur;

    Mass      const m_cutOffOxid       =
      m_ftEndOxid - ThrOxidRate    * SC::Stage3ThrBurnDur;

  private:
    // Consts used in MoI Optimisation:
    MoI             m_fuelMidMoIs[3];   // MoIs of Fuel in FuelTankMid
    Len             m_fuelMidCoM [3];   // CoM  ...

    MoI             m_fuelLowMoIs[3];   // MoIs of Fuel in FuelTankLow
    Len             m_fuelLowCoM [3];   // CoM  ...

    MoI             m_oxidMidMoIs[3];   // MoIs of Oxid in OxidTankMid
    Len             m_oxidMidCoM [3];   // CoM  ...

    MoI             m_oxidLowMoIs[3];   // MoIs of Oxid in OxidTankLow
    Len             m_oxidLowCoM [3];   // CoM  ...
 
  public:
    //=======================================================================//
    // Ctors:                                                                //
    //=======================================================================//
    // Default / Non-Default Ctor:
    Soyuz21b_Stage3();

    // Copy Ctor is auto-generated (for clarity):
    Soyuz21b_Stage3(Soyuz21b_Stage3 const&) = default;

    //=======================================================================//
    // Accessors:                                                            //
    //=======================================================================//
    TrCone    const& GetForeSection() const { return m_foreSection; }
    SpherSegm const& GetFuelTankUp () const { return m_fuelTankUp;  }
    TrCone    const& GetFuelTankMid() const { return m_fuelTankMid; }
    SpherSegm const& GetFuelTankLow() const { return m_fuelTankLow; }
    TrCone    const& GetEquipBay   () const { return m_equipBay;    }
    SpherSegm const& GetOxidTankUp () const { return m_oxidTankUp;  }
    TrCone    const& GetOxidTankMid() const { return m_oxidTankMid; }
    SpherSegm const& GetOxidTankLow() const { return m_oxidTankLow; }

    //=======================================================================//
    // Dynamic Params as functions of Flight Time:                           //
    //=======================================================================//
    StageDynParams GetDynParams(Time a_t) const;
  };
}
