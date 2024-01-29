// vim:ts=2:et
//===========================================================================//
//                            "Soyuz21b_Stage3.hpp":                         //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Soyuz21b_Stage3.h"
#include "SpaceBallistics/Soyuz21b_Consts.h"
#include "SpaceBallistics/Soyuz21b_Head_Consts.h"
#include "SpaceBallistics/Soyuz21b_Stage3_Consts.h"
#include "SpaceBallistics/Propellants.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace S3C = Soyuz21b_Stage3_Consts;
  namespace SHC = Soyuz21b_Head_Consts;
  namespace SC  = Soyuz21b_Consts;

  //=========================================================================//
  // Default / Non-Default Ctor:                                             //
  //=========================================================================//
  Soyuz21b_Stage3::Soyuz21b_Stage3()
  : //-----------------------------------------------------------------------//
    // ForeSection: Cylinder, UnKnownMass                                    //
    //-----------------------------------------------------------------------//
    m_foreSection
    (
      S3C::ForeX0,
      S3C::D,
      S3C::ForeH,
      Density(0.0)                  // No Propellant in this section
    ),

    //-----------------------------------------------------------------------//
    // Fuel Tank:                                                            //
    //-----------------------------------------------------------------------//
    // Up: HemiSphere, UnKnownMass
    m_fuelTankUp
    (
      false,                        // Facing Up = Left = !Right
      m_foreSection.GetRight()[0],  // Bottom of ForeSection
      S3C::D,
      RG1Dens                       // Naftil
    ),
    // Mid: Cylinder, UnKnownMass
    m_fuelTankMid
    (
      m_foreSection.GetRight()[0],  // Common base with FuelTankUp
      S3C::D,
      S3C::FuelTankMidH,
      RG1Dens                       // Naftil
    ),
    // Low: SpherSegm, UnKnownMass
    m_fuelTankLow
    (
      true,                         // Facing Down = Right
      m_fuelTankMid.GetRight()[0],  // Common base with FuelTankMid (right)
      S3C::D,
      S3C::FuelTankLowH,
      RG1Dens                       // Naftil
    ),

    //-----------------------------------------------------------------------//
    // Equipment Bay (containing the Control System): Cylinder, UnKnownMass  //
    //-----------------------------------------------------------------------//
    m_equipBay
    (
      m_fuelTankMid.GetRight()[0],  // Common base with FuelTankMid (right)
      S3C::D,
      S3C::EquipBayH,
      Density(0.0)                  // No Propellant in this Section
    ),

    //-----------------------------------------------------------------------//
    // Oxidiser Tank:                                                        //
    //-----------------------------------------------------------------------//
    // Up: HemiSpehere, UnKnownMass
    m_oxidTankUp
    (
      false,                        // Facing Up = Left = !Right
      m_equipBay.GetRight()[0],     // Common base with EquipBay (right)
      S3C::D,
      LOxDens
    ),
    // Mid: Cylinder, UnKnownMass
    m_oxidTankMid
    (
      m_equipBay.GetRight()[0],     // Common base with EquipBay (right)
      S3C::D,
      S3C::OxidTankMidH,
      LOxDens
    ),
    // Low: HemiSpehere, UnKnownMass
    m_oxidTankLow
    (
      true,                         // Facing Down = Right
      m_oxidTankMid.GetRight()[0],  // Common base with OxidTankMid (right)
      S3C::D,
      LOxDens
    ),

    //-----------------------------------------------------------------------//
    // Aft Section (jettisonable): Cylinder                                  //
    //-----------------------------------------------------------------------//
    m_aftSection
    (
      m_oxidTankMid.GetRight()[0],  // Common vase with OxidTankMid (right)
      S3C::D,
      S3C::AftH,
      Density(0.0),                 // No Propellant in this Section
      S3C::AftMass
    ),

    //-----------------------------------------------------------------------//
    // RD-0124 Engine, modeled as a PointMass:                               //
    //-----------------------------------------------------------------------//
    m_engine
    (
      m_oxidTankMid.GetRight()[0] + S3C::EngCoMdX,
      0.0_m,
      0.0_m,
      S3C::EngMass
    ),

    //-----------------------------------------------------------------------//
    // For Information Only:                                                 //
    //-----------------------------------------------------------------------//
    m_maxFuelMass
    (
      m_fuelTankUp .GetPropCap() + m_fuelTankMid.GetPropCap() +
      m_fuelTankLow.GetPropCap()
    ),

    m_maxOxidMass
    (
      m_oxidTankUp .GetPropCap() + m_oxidTankMid.GetPropCap() +
      m_oxidTankLow.GetPropCap()
    ),

    // NB: The following takes into account 2 periods of Throttled run as well:
    m_maxBurnTimeFT
    (
      std::min((S3C::FuelMass - S3C::FuelRem) / S3C::FuelRateFT,
               (S3C::OxidMass - S3C::OxidRem) / S3C::OxidRateFT)
      - 2.0 * SC::Stage3ThrBurnDur * S3C::Throttling
    )
  {
    //-----------------------------------------------------------------------//
    // Geometry Checks:                                                      //
    //-----------------------------------------------------------------------//
    // The bottom of the FuelTank must not touch the top of the Oxidiser Tank:
    assert(m_fuelTankLow.GetRight()[0] < m_oxidTankUp.GetLeft()[0]);

    // On the other hand, the top of the FuelTank appears over the over-all
    // cylindrical shell (but inside both kinds of InterStages):
    Len fuelTankTop = m_fuelTankUp.GetLeft()[0];
    assert(IsNeg(fuelTankTop) &&
           fuelTankTop > - SHC::InterStageLargeH &&
           fuelTankTop > - SHC::InterStageSmallH);
           
    // Also, the bottom of the OxidiserTank must be well inside the over-all
    // cylindrical shell:
    assert(m_oxidTankLow.GetRight()[0] < S3C::H);

    //-----------------------------------------------------------------------//
    // Mass and Time Checks:                                                 //
    //-----------------------------------------------------------------------//
    assert(S3C::FuelMass < m_maxFuelMass && S3C::OxidMass < m_maxOxidMass);

    // Max burn time at full thrust:
    assert(SC::Stage3ThrottlTime - SC::Stage3FullThrustTime
           < m_maxBurnTimeFT);

    // Needless to say, AftJet must occur before FullThrustEnd:
    assert(SC::Stage3AftJetTime < SC::Stage3FullThrustTime +
                                  SC::Stage3FTBurnDur);

    //-----------------------------------------------------------------------//
    // Set up Masses for all CEs:                                            //
    //-----------------------------------------------------------------------//
    ConstrElement::SetTotalMass
    (
      S3C::EmptyMass - S3C::AftMass - S3C::EngMass,

      { &m_foreSection,
        &m_fuelTankUp,  &m_fuelTankMid, &m_fuelTankLow,
        &m_equipBay,
        &m_oxidTankUp,  &m_oxidTankMid, &m_oxidTankLow
      }
    );

    //-----------------------------------------------------------------------//
    // Compute MoIs and CoMs of Fuel and Oxid sections:                      //
    //-----------------------------------------------------------------------//
    // NB: MoIs and CoMs are for Liquids only, w/o Tank walls (WithEmpty=false):
    //
/*
    m_fuelTankMid(m_fuelTankMid.GetPropCap(), false,
                  m_fuelMidMoIs, m_fuelMidCoM);
    m_fuelTankLow(m_fuelTankLow.GetPropCap(), false,
                  m_fuelLowMoIs, m_fuelLowCoM);

    m_oxidTankMid(m_oxidTankMid.GetPropCap(), false,
                  m_oxidMidMoIs, m_oxidMidCoM);
    m_oxidTankLow(m_oxidTankLow.GetPropCap(), false,
                  m_oxidLowMoIs, m_oxidLowCoM);
*/
  }

  //=========================================================================//
  // "GetDynParams":                                                         //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams Soyuz21b_Stage3::GetDynParams(Time a_t) const
  {
    StageDynParams res;           // XXX: Uninitialised yet

    //-----------------------------------------------------------------------//
    // Masses:                                                               //
    //-----------------------------------------------------------------------//
    Mass fullMass(NaN<double>);
    Mass fuelMass(NaN<double>);
    Mass oxidMass(NaN<double>);

    if (a_t < IgnTime)        // Before Stage3 Ignition
    {
      fullMass = FullMass;
      fuelMass = S3C::FuelMass;
      oxidMass = S3C::OxidMass;
    }
    else
    if (a_t < FTStartTime)    // Running at Initially-Throttled Thrust
    {
      Time   dt = a_t - IgnTime;
      fullMass  = FullMass      - ThrMassRate * dt;
      fuelMass  = S3C::FuelMass - ThrFuelRate * dt;
      oxidMass  = S3C::OxidMass - ThrOxidRate * dt;
    }
    else
    if (a_t < AftJetTime)     // Running at Full Thrust, with Aft (yet)
    {
      Time   dt = a_t - FTStartTime;
      fullMass  = FTStartMass - S3C::MassRateFT * dt;
      fuelMass  = FTStartFuel - S3C::FuelRateFT * dt;
      oxidMass  = FTStartOxid - S3C::OxidRateFT * dt;
    }
    else
    if (a_t < m_ftEndTime)    // Running at Full Thrust, no Aft
    {
      Time  dt  = a_t - AftJetTime;
      Time  dt1 = a_t - FTStartTime;
      fullMass  = AftJetMass  - S3C::MassRateFT * dt;
      fuelMass  = FTStartFuel - S3C::FuelRateFT * dt1;
      oxidMass  = FTStartOxid - S3C::OxidRateFT * dt1;
    }
    else
    if (a_t < m_cutOffTime)   // Running at Throttled Thrust again
    {
      Time dt   = a_t - m_ftEndTime;
      fullMass  = m_ftEndMass - ThrMassRate * dt;
      fuelMass  = m_ftEndFuel - ThrFuelRate * dt;
      oxidMass  = m_ftEndOxid - ThrOxidRate * dt;
    }
    else                      // After the Engine Cut-Off: "Spent" Stage
    {
      fullMass  = m_cutOffMass;
      fuelMass  = m_cutOffFuel;
      oxidMass  = m_cutOffOxid;
    }

    // Save the masses in the "res":
    res.m_fullMass = fullMass;
    res.m_fuelMass = fuelMass;
    res.m_oxidMass = oxidMass;

/*
    //-----------------------------------------------------------------------//
    // Moments of Inertia and Center of Masses:                              //
    //-----------------------------------------------------------------------//
    // XXX: Gases used for Stage3 tanks pressurisation are NOT included yet
    // into the MoI computation:
    //-----------------------------------------------------------------------//
    // Fuel:                                                                 //
    //-----------------------------------------------------------------------//
    if (fuelMass > m_fuelTankLow.GetPropCap() + m_fuelTankMod.GetPropCap())
    {
      // Fuel level is within the FuelTankUp:
      Mass fuelMassUp =
        fuelMass - m_fuelTankLow.GetPropCap() - m_fuelTankMod.GetPropCap();

      MoI fuelUpMoIs[3];
      Len fuelUpCoM [3];

      m_fuelTankUp(fuelMassUp, false, fuelUpMoIs, fuelUpCoM);

      // Add the MoIs and CoM from the Mid and Low sections. For that, create
      // a tmp "ConstrElement":
      ConstrElement tmp(Area(0.0), 

        res.m_mois[i] += m_fuelMidMoIs[i] + m_fuelLowMoIs[i];

      }

    }
*/
    return res;
  }
}
