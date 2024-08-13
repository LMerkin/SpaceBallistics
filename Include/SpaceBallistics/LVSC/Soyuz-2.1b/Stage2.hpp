// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.hpp":               //
//         Mathematical Model of the "Soyuz-2.1b" Stage2 ("Block A")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#include "SpaceBallistics/Utils.hpp"

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
    Time                   a_t,
    Pressure               a_p,      // Curr Atmospheric Pressure
    VernDeflections const& a_vern_defls
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    // We currently do not allow any times prior to LiftOff:
    assert(!IsNeg(a_t));

    // Also check the ShutDown sequnce: Verniers thtottle first, then the Main
    // engine:
    static_assert(VernThrottlTime < MainThrottlTime);

    // The atmospheric pressure is >= 0 obviously:
    assert(!IsNeg(a_p));

    // Check the MassRates at FullThrust:
    static_assert((OxidMR + FuelMR + H2O2MR).ApproxEquals(FullMR));

    StageDynParams<LVSC::Soyuz21b> res;               // XXX: All-0s yet...

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    // XXX: We assume for simplicity that both "Vac" and "SL" Thrust values are
    // throttled in the same proportion:

    Mass      fuelMass         (NaN<double>);
    Mass      oxidMass         (NaN<double>);
    Mass      h2o2Mass         (NaN<double>);
    MassRate  fuelMassDot      (NaN<double>);
    MassRate  oxidMassDot      (NaN<double>);
    MassRate  h2o2MassDot      (NaN<double>);

    Force     absMainThrustVac (NaN<double>);         // Main Engine
    Force     absMainThrustSL  (NaN<double>);
    Force     absVernThrustVac1(NaN<double>);         // Each Vernier Chamber
    Force     absVernThrustSL1 (NaN<double>);


    if (a_t < VernThrottlTime)
    {
      // Full-Thrust Mode:
      fuelMass          = FuelMass0 - FuelMR  * a_t;
      oxidMass          = OxidMass0 - OxidMR  * a_t;
      h2o2Mass          = H2O2Mass0 - H2O2MR  * a_t;
      fuelMassDot       = - FuelMR;
      oxidMassDot       = - OxidMR;
      h2o2MassDot       = - H2O2MR;

      absMainThrustVac  = ThrustMainVac;
      absMainThrustSL   = ThrustMainSL;
      absVernThrustVac1 = ThrustVernVac1;
      absVernThrustSL1  = ThrustVernSL1;
    }
    else
    if (a_t < MainThrottlTime)
    {
      // Verniers are throttled, but the Main Engine is not yet:
      Time  dt          = a_t - VernThrottlTime;
      fuelMass          = FuelMassV - FuelMRV * dt;
      oxidMass          = OxidMassV - OxidMRV * dt;
      h2o2Mass          = H2O2Mass0 - H2O2MR  * a_t;
      fuelMassDot       = - FuelMRV;
      oxidMassDot       = - OxidMRV;
      h2o2MassDot       = - H2O2MR;

      absMainThrustVac  = ThrustMainVac;
      absMainThrustSL   = ThrustMainSL;
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
    else
    if (a_t < CutOffTime)
    {
      Time  dt          = a_t - MainThrottlTime;
      fuelMass          = FuelMassM - FuelMRM * dt;
      oxidMass          = OxidMassM - OxidMRM * dt;
      h2o2Mass          = H2O2Mass0 - H2O2MR  * a_t;
      fuelMassDot       = - FuelMRM;
      oxidMassDot       = - OxidMRM;
      h2o2MassDot       = - H2O2MR;

      absMainThrustVac  = ThrustMainVac  * ShutDownThrottlLevel;
      absMainThrustSL   = ThrustMainSL   * ShutDownThrottlLevel;
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
    else
    {
      fuelMass          = FuelMassC;
      oxidMass          = OxidMassC;
      h2o2Mass          = H2O2MassC;
      fuelMassDot       = MassRate(0.0);
      oxidMassDot       = MassRate(0.0);
      h2o2MassDot       = MassRate(0.0);

      absMainThrustVac  = Force   (0.0);
      absMainThrustSL   = Force   (0.0);
      absVernThrustVac1 = Force   (0.0);
      absVernThrustSL1  = Force   (0.0);
    }
    // "fullMass" is "FullMass0" (at lift-off) minus the mass of Fuel, Oxid and
    // H2O2 spent:
    Mass  fuelSpent     = FuelMass0 - fuelMass;
    Mass  oxidSpent     = OxidMass0 - oxidMass;
    Mass  h2o2Spent     = H2O2Mass0 - h2o2Mass;
    assert(!(IsNeg(fuelSpent) || IsNeg(oxidSpent) || IsNeg(h2o2Spent)));

    Mass  fullMass      = FullMass0 - (fuelSpent + oxidSpent + h2o2Spent);

    // Check the Curr Masses:
    assert(fullMass >= MinEndMass && fuelMass >= FuelRem &&
           oxidMass >= OxidRem    && h2o2Mass >= H2O2Rem);

    // If OK: Save the masses in the "res":
    res.m_fullMass  = fullMass;
    res.m_fuelMass  = fuelMass;
    res.m_oxidMass  = oxidMass;

    //-----------------------------------------------------------------------//
    // Thrust Vector:                                                        //
    //-----------------------------------------------------------------------//
    // Adjust the Thrust values to the curr atmospheric pressure:
    //
    double slC  = double(a_p / p0);    // 1 when we are @ SL
    assert(0.0 <= slC && slC <= 1.0);
    double vacC = 1.0 - slC;           // 1 when we are in Vac

    Force absMainThrust  = vacC * absMainThrustVac  + slC * absMainThrustSL;
    Force absVernThrust1 = vacC * absVernThrustVac1 + slC * absVernThrustSL1;

    Force thrustX = absMainThrust;
    Force thrustY;  // Initially 0
    Force thrustZ;  // ditto

    // Consider VernGimbalAngles of all 4 Vernier Chambers.  The formulas here
    // are the same as for Stage3, because the Stage2 Verniers and Stage3 Main
    // (Gimbaled) Chambers are installed in the same planes:
    //
    for (size_t i = 0; i < 4; ++i)
    {
      Angle_deg  A   = a_vern_defls[i];
      assert(Abs(A) <= VernGimbalAmpl);

      double sinA = Sin(double(To_Angle(A)));
      double cosA = Cos(double(To_Angle(A)));
      thrustX    += absVernThrust1 * cosA;

      // NB: The ThrustVector rotation in the YZ plane is OPPOSITE to the
      // corresp Vernier Deflection:
      switch (i)
      {
      case 0:
        thrustZ += absVernThrust1 * sinA;
        break;
      case 1:
        thrustY -= absVernThrust1 * sinA;
        break;
      case 2:
        thrustZ -= absVernThrust1 * sinA;
        break;
      case 3:
        thrustY += absVernThrust1 * sinA;
        break;
      default:
        assert(false);
      }
    }
    res.m_thrust = ME::ForceVE{{thrustX, thrustY, thrustZ}};

/*
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
