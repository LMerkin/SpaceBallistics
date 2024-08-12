// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.hpp":               //
//         Mathematical Model of the "Soyuz-2.1b" Stage2 ("Block A")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#include <stdexcept>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage2::GetDynParams":                                        //
  //=========================================================================//
  // "a_t" is Flight Time since the "Contact Separation" event:
  //
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage2::GetDynParams
  (
    Time                      a_t,
    Pressure                  a_p,      // Curr Atmospheric Pressure
    ChamberDeflections const& a_chamber_defls
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    // We currently do not allow any times prior to LiftOff:
    if (UNLIKELY(IsNeg(a_t)))
      throw std::logic_error("Soyuz21b_Stage2: t < 0 not allowed");

    // Also check the ShutDown sequnce: Verniers thtottle first, then the Main
    // engine:
    static_assert(VernThrottlTime < MainThrottlTime);

    // The atmospheric pressure is >= 0 obviously:
    assert(!IsNeg(a_p));

    // Check the MassRates at FullThrust:
    static_assert((OxidMR + FuelMR + H2O2MR).ApproxEquals(TotalMR));

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    StageDynParams<LVSC::Soyuz21b> res;               // XXX: All-0s yet...

    Mass      oxidMass         (NaN<double>);
    Mass      fuelMass         (NaN<double>);
    Mass      h2o2Mass         (NaN<double>);
    Mass      fullMass         (NaN<double>);

    MassRate  oxidMassDot      (NaN<double>);
    MassRate  fuelMassDot      (NaN<double>);
    MassRate  h2o2MassDot      (NaN<double>);

    Force     absMainThrustVac (NaN<double>);         // Main Engine
    Force     absVernThrustVac1(NaN<double>);         // Each Vernier Chamber

    if (a_t < VernThrottlTime)
    {
      // Full-Thrust Mode:
      oxidMass          = OxidMass0 - OxidMR  * a_t;  // Not including H2O2MR
      fuelMass          = FuelMass0 - FuelMR  * a_t;  //
      h2o2Mass          = H2O2Mass0 - H2O2MR  * a_t;
      fullMass          = FullMass0 - TotalMR * a_t;  // Including H2O2MR
      oxidMassDot       = - OxidMR;
      fuelMassDot       = - FuelMR;
      h2o2MassDot       = - H2O2MR;
//    absMainThrustVac  =
    }
    else
    if (a_t < MainThrottlTime)
    {
    }
    else
    if (a_t < CutOffTime)
    {
    }
    else
    {
    }

/*
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
              ("Soyuz21b_Stage2::GetDynParams: GimbalAmpl exceeded");

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
*/

    // All Done:
    return res;
  }
}
// End namespace SpaceBallistics
