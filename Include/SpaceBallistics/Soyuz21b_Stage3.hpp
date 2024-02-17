// vim:ts=2:et
//===========================================================================//
//                            "Soyuz21b_Stage3.hpp":                         //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Soyuz21b_Stage3.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "GetDynParams":                                                         //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams Soyuz21b_Stage3::GetDynParams(Time a_t)
  {
    StageDynParams res;           // XXX: Uninitialised yet

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    Mass  fullMass(NaN<double>);
    Mass  fuelMass(NaN<double>);
    Mass  oxidMass(NaN<double>);
    Force thrust  (NaN<double>);

    if (a_t < IgnTime)        // Before Stage3 Ignition
    {
      fullMass = FullMass;
      fuelMass = FuelMass;
      oxidMass = OxidMass;
      thrust   = Force(0.0);
    }
    else
    if (a_t < FTStartTime)    // Running at Initially-Throttled Thrust
    {
      Time   dt = a_t - IgnTime;
      fullMass  = FullMass  - ThrMassRate * dt;
      fuelMass  = FuelMass  - ThrFuelRate * dt;
      oxidMass  = OxidMass  - ThrOxidRate * dt;
      thrust    = ThrottlThrust;
    }
    else
    if (a_t < AftJetTime)     // Running at Full Thrust, with Aft (yet)
    {
      Time   dt = a_t - FTStartTime;
      fullMass  = FTStartFullMass - MassRateFT * dt;
      fuelMass  = FTStartFuelMass - FuelRateFT * dt;
      oxidMass  = FTStartOxidMass - OxidRateFT * dt;
      thrust    = ThrustVac;
    }
    else
    if (a_t < FTEndTime)      // Running at Full Thrust, no Aft
    {
      Time  dt1 = a_t - FTStartTime;
      Time  dt2 = a_t - AftJetTime;
      fullMass  = AftJetFullMass  - MassRateFT * dt2;
      fuelMass  = FTStartFuelMass - FuelRateFT * dt1;
      oxidMass  = FTStartOxidMass - OxidRateFT * dt1;
      thrust    = ThrustVac;
    }
    else
    if (a_t < CutOffTime)     // Running at Throttled Thrust again
    {
      Time dt   = a_t - FTEndTime;
      fullMass  = FTEndFullMass   - ThrMassRate * dt;
      fuelMass  = FTEndFuelMass   - ThrFuelRate * dt;
      oxidMass  = FTEndOxidMass   - ThrOxidRate * dt;
      thrust    = ThrottlThrust;
    }
    else                      // After the Engine Cut-Off: "Spent" Stage
    {
      fullMass  = CutOffFullMass;
      fuelMass  = CutOffFuelMass;
      oxidMass  = CutOffOxidMass;
      thrust    = Force(0.0);
    }

    // Verify the Masses:
    DEBUG_ONLY
    (
      //  "egMass" includes EmptyMass and GasesMass:
      Mass egMass = fullMass - fuelMass - oxidMass;
      assert
        (egMass.ApproxEquals(a_t < AftJetTime ? EGMassBefore : EGMassAfter));
    )
    // Save the masses in the "res":
    res.m_fullMass = fullMass;
    res.m_fuelMass = fuelMass;
    res.m_oxidMass = oxidMass;
    res.m_thrust   = thrust;

    //-----------------------------------------------------------------------//
    // Moments of Inertia and Center of Masses:                              //
    //-----------------------------------------------------------------------//
    // XXX: For the purpose of MoI computation, we include the masses  of Tank
    // Pressurisation Gases into the corresp Fuel/Oxid masses (proportional to
    // the volumes of the corresp Tanks):
    //
    // Fuel:
    assert(IsPos(fuelMass) && fuelMass <= FuelMass && FuelMass < FuelTankMC);
    CE fuelCE =
      (fuelMass > FuelTankLowMidMC)
      ? // Fuel level is within the FuelTankUp:
        FuelTankUp .GetPropCE(fuelMass - FuelTankLowMidMC) + FuelLowMidCE
      :
      (fuelMass > FuelTankLowMC)
      ? // Fuel level is within the FuelTankMid:
        FuelTankMid.GetPropCE(fuelMass - FuelTankLowMC)    + FuelLowCE
      :
        // Fuel level is within the FuelTankLow:
        FuelTankLow.GetPropCE(fuelMass);

    // Oxid:
    assert(IsPos(oxidMass) && oxidMass <= OxidMass && OxidMass < OxidTankMC);
    CE oxidCE =
      (oxidMass > OxidTankLowMidMC)
      ? // Oxid level is within the OxidTankUp:
        OxidTankUp .GetPropCE(oxidMass - OxidTankLowMidMC) + OxidLowMidCE
      :
      (oxidMass > OxidTankLowMC)
      ? // Oxid level is within the OxidTankMid:
        OxidTankMid.GetPropCE(oxidMass - OxidTankLowMC)    + OxidLowCE
      :
        // Oxid level is within the OxidTankLow:
        OxidTankLow.GetPropCE(oxidMass);

    // Full = (Empty+Gases) + Fuel + Oxid:
    CE fullCE = (a_t < AftJetTime ? EGBeforeCE : EGAfterCE) + fuelCE + oxidCE;

    // Double-check the Masses:
    assert(fullMass.ApproxEquals(fullCE.GetMass()) &&
           fuelMass.ApproxEquals(fuelCE.GetMass()) &&
           oxidMass.ApproxEquals(oxidCE.GetMass()));

    // Extract the the CoM and the MoIs:
    CE::Point const& com  = fullCE.GetCoM ();
    CE::MoIs  const& mois = fullCE.GetMoIs();

    res.m_com [0] =  com [0];
    res.m_com [1] =  com [1];
    res.m_com [2] =  com [2];

    res.m_mois[0] =  mois[0];
    res.m_mois[1] =  mois[1];
    res.m_mois[2] =  mois[2];

    // All Done:
    return res;
  }
}
