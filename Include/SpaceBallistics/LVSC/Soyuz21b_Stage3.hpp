// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/LVSC/Soyuz21b_Stage3.hpp":                //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/Soyuz21b_Stage3.h"
#include "SpaceBallistics/Utils.hpp"
#include <stdexcept>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage3::GetDynParams":                                        //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage3::GetDynParams(Time a_t, ThrustVecCtl const& a_thrust_ctl)
  {
    StageDynParams<LVSC::Soyuz21b> res;   // XXX: Empty (all 0s) yet...

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

    //-----------------------------------------------------------------------//
    // Thrust Vector:                                                        //
    //-----------------------------------------------------------------------//
    Force thrustX;        // Initially 0
    Force thrustY;        //
    Force thrustZ;        //
    Force chamberThrust = absThrust / 4.0;

    // Consider GimbalAngles of all 4 Chambers:
    for (size_t i = 0; i < 4; ++i)
    {
      Angle_deg A = a_thrust_ctl[i];
      if (Abs(A) > GimbalAmpl)
        throw std::invalid_argument
              ("Soyuz21b_Stage3::GetDynParams: GimbalAmpl exceeded");

      double sinA = Sin(double(To_Angle(A)));
      double cosA = Cos(double(To_Angle(A)));
      thrustX    += chamberThrust * cosA;

      // NB: The ThrustVector rotation in the YZ plane is opposite to the
      // corresp Chamber (more precisely, Nozzle) rotation:
      switch (i)
      {
      case 0:
        thrustZ -= chamberThrust * sinA;
        break;
      case 1:
        thrustY += chamberThrust * sinA;
        break;
      case 2:
        thrustZ += chamberThrust * sinA;
        break;
      case 3:
        thrustY -= chamberThrust * sinA;
        break;
      default:
        assert(false);
      }
    }
    res.m_thrust = ME::ForceVE{{thrustX, thrustY, thrustZ}};

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
    res.m_com  = fullME.GetCoM ();
    res.m_mois = fullME.GetMoIs();

    // All Done:
    return res;
  }
}
// End namespace SpaceBallistics
