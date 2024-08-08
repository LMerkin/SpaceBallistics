// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.hpp":               //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Utils.hpp"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"
#include <stdexcept>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage3::GetDynParams":                                        //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage3::GetDynParams
  (
    Time                      a_t,
    ChamberDeflections const& a_chamber_defls
  )
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
    if (a_t < CutOffTime)     // Running at Full Thrust, w/o Aft
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

      // Also, the Fuel and Oxid masses must not be below the physical low
      // limits:
      if (UNLIKELY(fuelMass < FuelRem || oxidMass < OxidRem))
        throw std::logic_error(std::format
             ("Soyuz21b_Stage3: t={}: Fuel and/or Oxid Mass too low: "
              "Fuel={} (Min={}), Oxid={} (Min={})",
              a_t, fuelMass, FuelRem, oxidMass, OxidRem));
    )
    // If OK: Save the masses in the "res":
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
      Angle_deg A = a_chamber_defls[i];
      if (Abs(A)  > GimbalAmpl)
        throw std::invalid_argument
              ("Soyuz21b_Stage3::GetDynParams: GimbalAmpl exceeded");

      double sinA = Sin(double(To_Angle(A)));
      double cosA = Cos(double(To_Angle(A)));
      thrustX    += chamberThrust * cosA;

      // NB: The ThrustVector rotation in the YZ plane is OPPOSITE to the
      // corresp Chamber Deflection:
      switch (i)
      {
      case 0:
        thrustZ += chamberThrust * sinA;
        break;
      case 1:
        thrustY -= chamberThrust * sinA;
        break;
      case 2:
        thrustZ -= chamberThrust * sinA;
        break;
      case 3:
        thrustY += chamberThrust * sinA;
        break;
      default:
        assert(false);
      }
    }
    res.m_thrust = ME::ForceVE{{thrustX, thrustY, thrustZ}};

    //-----------------------------------------------------------------------//
    // Moments of Inertia and Center of Masses:                              //
    //-----------------------------------------------------------------------//
    // NB:
    // (*) "fuelLevel" and "oxidLevel" are hooks provided for debugging purposes
    //     only;
    // (*) Gases (primarily He) are considered along with the EmptyMass; redist-
    //     ribution of Gases within Stage3 is NOT taken into account, whereas in
    //     reality, He is used for Fuel and Oxid tanks pressurisation, as fills
    //     the "empty" space above the remaining Fuel and Oxid.   However, this
    //     effect is considered negligible for Stage3 (though may have some imp-
    //     act on the Moments of Inertia for Stage2):
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
