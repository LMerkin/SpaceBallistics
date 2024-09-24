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
    Mass      liqN2Mass        (NaN<double>);

    MassRate  fuelMassDot      (NaN<double>);
    MassRate  oxidMassDot      (NaN<double>);
    MassRate  h2o2MassDot      (NaN<double>);
    MassRate  liqN2MassDot     (NaN<double>);

    Force     absMainThrustVac (NaN<double>);         // Main Engine
    Force     absMainThrustSL  (NaN<double>);
    Force     absVernThrustVac1(NaN<double>);         // Each Vernier Chamber
    Force     absVernThrustSL1 (NaN<double>);

    if (a_t < ThrottlTime)
    {
      // Full-Thrust Mode:
      fuelMass          = FuelMass0  - FuelMR  * a_t;
      oxidMass          = OxidMass0  - OxidMR  * a_t;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * a_t;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * a_t;
      fuelMassDot       = - FuelMR;
      oxidMassDot       = - OxidMR;
      h2o2MassDot       = - H2O2MR;
      liqN2MassDot      = - LiqN2MR;

      absMainThrustVac  = ThrustMainVac;
      absMainThrustSL   = ThrustMainSL;
      absVernThrustVac1 = ThrustVernVac1;
      absVernThrustSL1  = ThrustVernSL1;
    }
    else
    if (a_t < MainCutOffTime)
    {
      // Throttled Mode (for both Main Engine and Verniers):
      Time  dt          = a_t - ThrottlTime;
      fuelMass          = FuelMassT  - FuelMRT * dt;
      oxidMass          = OxidMassT  - OxidMRT * dt;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * a_t;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * a_t;
      fuelMassDot       = - FuelMRT;
      oxidMassDot       = - OxidMRT;
      h2o2MassDot       = - H2O2MR;
      liqN2MassDot      = - LiqN2MR;

      absMainThrustVac  = ThrustMainVac  * ShutDownThrottlLevel;
      absMainThrustSL   = ThrustMainSL   * ShutDownThrottlLevel;
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
    else
    if (a_t < CutOffTime)
    {
      // Only the Verniers continue operation, in the Throttled Mode:
      Time  dt          = a_t - MainCutOffTime;
      fuelMass          = FuelMassM  - FuelMRM * dt;
      oxidMass          = OxidMassM  - OxidMRM * dt;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * a_t;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * a_t;
      fuelMassDot       = - FuelMRM;
      oxidMassDot       = - OxidMRM;
      h2o2MassDot       = - H2O2MR;
      liqN2MassDot      = - LiqN2MR;

      absMainThrustVac  = Force   (0.0);
      absMainThrustSL   = Force   (0.0);
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
    else
    {
      fuelMass          = FuelMassC;
      oxidMass          = OxidMassC;
      h2o2Mass          = H2O2MassC;
      liqN2Mass         = LiqN2MassC;
      fuelMassDot       = MassRate(0.0);
      oxidMassDot       = MassRate(0.0);
      h2o2MassDot       = MassRate(0.0);
      liqN2MassDot      = MassRate(0.0);

      absMainThrustVac  = Force   (0.0);
      absMainThrustSL   = Force   (0.0);
      absVernThrustVac1 = Force   (0.0);
      absVernThrustSL1  = Force   (0.0);
    }
    //
    // "fullMass" is "FullMass0" (at LiftOff) minus the mass of Fuel, Oxid and
    // H2O2 spent (Once again, N2O2 is NOT spent):
    Mass  fuelSpent     = FuelMass0 - fuelMass;
    Mass  oxidSpent     = OxidMass0 - oxidMass;
    Mass  h2o2Spent     = H2O2Mass0 - h2o2Mass;
    assert(!(IsNeg(fuelSpent) || IsNeg(oxidSpent) || IsNeg(h2o2Spent)));

    Mass  fullMass      = FullMass0 - (fuelSpent + oxidSpent + h2o2Spent);

    // Check the Curr Masses and MassRates:
    assert(MinEndMass <= fullMass  && fullMass  <= FullMass0  &&
           FuelRem    <= fuelMass  && fuelMass  <= FuelMass0  &&
           OxidRem    <= oxidMass  && oxidMass  <= OxidMass0  &&
           H2O2Rem    <= h2o2Mass  && h2o2Mass  <= H2O2Mass0  &&
           LiqN2Rem   <= liqN2Mass && liqN2Mass <= LiqN2Mass0 &&
           !(IsPos(fuelMassDot) || IsPos(oxidMassDot) || IsPos(h2o2MassDot) ||
             IsPos(liqN2MassDot)));

    // If OK: Save the masses in the "res"; XXX: "h2o2Mass" is not installed --
    // there is no field in "res" for it; it is only used in the MoIs and CoM
    // computations below:
    res.m_fullMass  = fullMass;
    res.m_fuelMass  = fuelMass;
    res.m_oxidMass  = oxidMass;

    //-----------------------------------------------------------------------//
    // Thrust Vector:                                                        //
    //-----------------------------------------------------------------------//
    // Adjust the Thrust values to the curr atmospheric pressure:
    //
    double slC  = double(a_p / p0);    // 1 when we are @ SL
    double vacC = 1.0 - slC;           // 1 when we are in Vac
    assert(0.0 <= slC && slC <= 1.0 && 0 <= vacC && vacC <= 1.0);

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

    //-----------------------------------------------------------------------//
    // Moments of Inertia and Center of Masses:                              //
    //-----------------------------------------------------------------------//
    // XXX: Similar to Stage3, "fuelLevel", "oxidLevel" and "h2o2Level" are
    // provided for debugging purposes only:
    //
    // Fuel:
    Len fuelLevel(NaN<double>);
    ME fuelME =
      (fuelMass > FuelTankBtmMidMC)
      ? // Fuel level is within the FuelTankTop:
        FuelTankTop.GetPropBulkME
          (fuelMass - FuelTankBtmMidMC,     fuelMassDot,  &fuelLevel) +
        FuelBtmMidME
      :
      (fuelMass > FuelTankBtmMC)
      ? // Fuel level is within the FuelTankMid:
        FuelTankMid.GetPropBulkME
          (fuelMass - FuelTankBtmMC,        fuelMassDot,  &fuelLevel) +
        FuelBtmME
      :
        // Fuel level is within the FuelTankBtm:
        FuelTankBtm.GetPropBulkME(fuelMass, fuelMassDot,  &fuelLevel);
    assert(IsPos(fuelLevel));

    // Oxid:
    Len oxidLevel(NaN<double>);

    ME oxidME =
      (oxidMass > OxidTankBtmLowUpMC)
      ? // Oxid level is within the OxidTankTop:
        OxidTankTop.GetPropBulkME
          (oxidMass - OxidTankBtmLowUpMC,   oxidMassDot,  &oxidLevel) +
        OxidBtmLowUpME
      :
      (oxidMass > OxidTankBtmLowMC)
      ? // Oxid level is within the OxidTankUp:
        OxidTankUp .GetPropBulkME
          (oxidMass - OxidTankBtmLowMC,     oxidMassDot,  &oxidLevel) +
        OxidBtmLowME
      :
      (oxidMass > OxidTankBtmMC)
      ? // Oxid level is within the OxidTankLow:
        OxidTankLow.GetPropBulkME
          (oxidMass - OxidTankBtmMC,        oxidMassDot,  &oxidLevel) +
        OxidBtmME
      :
        // Oxid level is within the OxidTankBtm:
        OxidTankBtm.GetPropBulkME(oxidMass, oxidMassDot,  &oxidLevel);
    assert(IsPos(oxidLevel));

    // H2O2:
    Len h2o2Level(NaN<double>);

    ME h2o2ME =
      (h2o2Mass > H2O2TankBtmMidMC)
      ? // H2O2 level is within the H2O2TankTop:
        H2O2TankTop.GetPropBulkME
          (h2o2Mass - H2O2TankBtmMidMC,     h2o2MassDot,  &h2o2Level) +
        H2O2BtmMidME
      :
      (h2o2Mass > H2O2TankBtmMC)
      ? // H2O2 level is within the H2H2TankMid:
        H2O2TankMid.GetPropBulkME
          (h2o2Mass - H2O2TankBtmMC,        h2o2MassDot,  &h2o2Level) +
        H2O2BtmME
      :
        // H2O2 level is within the H2O2TankBtm:
        H2O2TankBtm.GetPropBulkME(h2o2Mass, h2o2MassDot,  &h2o2Level);
    assert(IsPos(h2o2Level));

    // LiqN2:
    Len liqN2Level(NaN<double>);

    ME liqN2ME =
      (liqN2Mass > LiqN2TankBtmMC)
      ? // LiqN2 level is within the LiqN2TankTop:
        LiqN2TankTop.GetPropBulkME
          (liqN2Mass - LiqN2TankBtmMC,        liqN2MassDot, &liqN2Level) +
        LiqN2BtmME
      :
        // LiqN2 level is within the LiqN2TankBtm:
        LiqN2TankBtm.GetPropBulkME(liqN2Mass, liqN2MassDot, &liqN2Level);
    assert(IsPos(liqN2Level));

    // GasN2:
    // Unlike Stage3, here we cannot neglect the effects  of the Tank Pressuri-
    // sation Gases (in this case, N2), because their mass is reasonable large.
    // N2 is NOT spent (exhausted), but its distribution within the Tanks chan-
    // ges over the time. Gasesous N2 is concentrated in the upper parts of the
    // Fuel and Oxid Tanks, gradually increasing in volume.
    // The "complementary" ("CME") parts of the Fuel and Oxid bulks (which are
    // actually filled with N2):
    //
    ME   fuelCME   = FuelME  - fuelME;
    ME   oxidCME   = OxidME  - oxidME;
    ME   complME   = fuelCME + oxidCME;

    // Mass of the Gaseous N2:
    Mass gasN2Mass = N2Mass - liqN2Mass;
    assert(GasN2Mass0 <= gasN2Mass && gasN2Mass < N2Mass - LiqN2Rem);
    ME   gasN2ME   =
      ME::ProRateMass<true>(complME, double(gasN2Mass / complME.GetMass()));
    /*
    // Density of the Gaseous N2 (enable for debugging only):
    Vol     gasN2Vol = complME.GetEnclVol();
    Density gasN2Rho = gasN2Mass / gasN2Vol;
    */
    // Full = Empty + N2 + Fuel + Oxid + H2O2 + LiqN2 + GasN2:
    ME fullME = EmptyME  + fuelME + oxidME + h2o2ME + liqN2ME + gasN2ME;

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
