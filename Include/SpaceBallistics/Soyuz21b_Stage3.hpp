// vim:ts=2:et
//===========================================================================//
//                            "Soyuz21b_Stage3.hpp":                         //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#include "SpaceBallistics/Soyuz21b_Stage3.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "GetDynParams":                                                         //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
/*
  StageDynParams Soyuz21b_Stage3::GetDynParams(Time a_t)
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
    return res;
  }
*/
}
