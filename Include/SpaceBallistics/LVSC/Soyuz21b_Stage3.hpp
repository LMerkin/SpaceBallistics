// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/LVSC/Soyuz21b_Stage3.hpp":                //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/Soyuz21b_Stage3.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage3::GetDynParams":                                        //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams Soyuz21b_Stage3::GetDynParams(Time a_t)
  {
    StageDynParams res;   // XXX: Uninitialised yet

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    using MR = ::SpaceBallistics::MassRate;

    Mass    fullMass   (NaN<double>);
    Mass    fuelMass   (NaN<double>);
    Mass    oxidMass   (NaN<double>);
    MR      fuelMassDot(NaN<double>);
    MR      oxidMassDot(NaN<double>);
    Force   absThrust  (NaN<double>);

    if (a_t < IgnTime)        // Before Stage3 Ignition
    {
      fullMass    = FullMass;
      fuelMass    = FuelMass;
      oxidMass    = OxidMass;
      absThrust   = Force(0.0);
      fuelMassDot = MR(0.0);
      oxidMassDot = MR(0.0);
    }
    else
    if (a_t < AftJetTime)     // Running at Full Thrust, with Aft (yet)
    {
      Time     dt = a_t      - IgnTime;
      fullMass    = FullMass - MassRate * dt;
      fuelMass    = FuelMass - FuelRate * dt;
      oxidMass    = OxidMass - OxidRate * dt;
      absThrust   = ThrustVac;
      fuelMassDot = - FuelRate;
      oxidMassDot = - OxidRate;
    }
    else
    if (a_t < CutOffTime)     // Running at Throttled Thrust again
    {
      Time dt0    = a_t             - IgnTime;
      Time dt1    = a_t             - AftJetTime;
      fullMass    = AftJetFullMass  - MassRate * dt1;
      fuelMass    = FuelMass        - FuelRate * dt0;
      oxidMass    = OxidMass        - OxidRate * dt0;
      absThrust   = ThrustVac;
      fuelMassDot = - FuelRate;
      oxidMassDot = - OxidRate;
    }
    else                      // After the Engine Cut-Off: "Spent" Stage
    {
      fullMass    = CutOffFullMass;
      fuelMass    = CutOffFuelMass;
      oxidMass    = CutOffOxidMass;
      absThrust   = Force(0.0);
      fuelMassDot = MR(0.0);
      oxidMassDot = MR(0.0);
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
    res.m_fullMass  = fullMass;
    res.m_fuelMass  = fuelMass;
    res.m_oxidMass  = oxidMass;
    // FIXME: We currently assume that the Thrust vector is parallel to the
    // embedded OX axis (pointing in the positive direction):
    res.m_thrust[0] = absThrust;
    res.m_thrust[1] = Force(0.0);
    res.m_thrust[2] = Force(0.0);

    //-----------------------------------------------------------------------//
    // Moments of Inertia and Center of Masses:                              //
    //-----------------------------------------------------------------------//
    // XXX: For the purpose of MoI computation, we include the masses  of Tank
    // Pressurisation Gases into the corresp Fuel/Oxid masses (proportional to
    // the volumes of the corresp Tanks):
    // XXX: "fuelLevel" and "oxidLevel" are hooks provided for debugging purpo-
    // ses only:
    //
    // Fuel:
    assert(IsPos(fuelMass) && fuelMass <= FuelMass && FuelMass < FuelTankMC);
    Len fuelLevel = 0.0_m;
    ME fuelME =
      (fuelMass > FuelTankLowMidMC)
      ? // Fuel level is within the FuelTankUp:
        FuelTankUp .GetPropBulkME
          (fuelMass - FuelTankLowMidMC, fuelMassDot, &fuelLevel) +
        FuelLowMidME
      :
      (fuelMass > FuelTankLowMC)
      ? // Fuel level is within the FuelTankMid:
        FuelTankMid.GetPropBulkME
          (fuelMass - FuelTankLowMC,    fuelMassDot, &fuelLevel) +
        FuelLowME
      :
        // Fuel level is within the FuelTankLow:
        FuelTankLow.GetPropBulkME(fuelMass, fuelMassDot, &fuelLevel);

    // Oxid:
    Len oxidLevel = 0.0_m;
    assert(IsPos(oxidMass) && oxidMass <= OxidMass && OxidMass < OxidTankMC);
    ME oxidME =
      (oxidMass > OxidTankLowMidMC)
      ? // Oxid level is within the OxidTankUp:
        OxidTankUp .GetPropBulkME
          (oxidMass - OxidTankLowMidMC, oxidMassDot, &oxidLevel) +
        OxidLowMidME
      :
      (oxidMass > OxidTankLowMC)
      ? // Oxid level is within the OxidTankMid:
        OxidTankMid.GetPropBulkME
          (oxidMass - OxidTankLowMC,    oxidMassDot, &oxidLevel) +
        OxidLowME
      :
        // Oxid level is within the OxidTankLow:
        OxidTankLow.GetPropBulkME(oxidMass, oxidMassDot, &oxidLevel);

    // Full = (Empty+Gases) + Fuel + Oxid:
    ME fullME = (a_t < AftJetTime ? EGBeforeME : EGAfterME) + fuelME + oxidME;

    // Double-check the Masses:
    assert(fullMass.ApproxEquals(fullME.GetMass()) &&
           fuelMass.ApproxEquals(fuelME.GetMass()) &&
           oxidMass.ApproxEquals(oxidME.GetMass()));

    // Extract the the CoM and the MoIs:
    ME::PosVE const& com  = fullME.GetCoM ();
    ME::MoITE const& mois = fullME.GetMoIs();

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
// End namespace SpaceBallistics
