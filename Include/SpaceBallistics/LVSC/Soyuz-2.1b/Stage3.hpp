// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.hpp":               //
//         Mathematical Model of the "Soyuz-2.1b" Stage3 ("Block I")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage3::GetDynParams":                                        //
  //=========================================================================//
  // "a_t" is the Flight Time since the LiftOff ("Contact Separation") event:
  //
  StageDynParams<LVSC::Soyuz21b> Soyuz21b_Stage3::GetDynParams
  (
    FlightTime                a_ft,
    ChamberDeflections const& a_chamber_defls
  )
  {
    // We currently do not allow any times prior to LiftOff:
    assert(SC::LiftOffTime <= a_ft);

    // The absolute curr time must be defined in this case:
    TT  tt = TT(a_ft);
    assert(!tt.IsUnDef());

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    Mass      fullMass   (NaN<double>);
    Mass      fuelMass   (NaN<double>);
    Mass      oxidMass   (NaN<double>);
    MassRate  fuelMassDot(NaN<double>);
    MassRate  oxidMassDot(NaN<double>);
    Force     absThrust  (NaN<double>);

    if (a_ft < IgnTime)      // Before Stage3 Ignition
    {
      fullMass    = FullMass;
      fuelMass    = FuelMass;
      oxidMass    = OxidMass;
      absThrust   = Force   (0.0);
      fuelMassDot = MassRate(0.0);
      oxidMassDot = MassRate(0.0);
    }
    else
    if (a_ft < AftJetTime)   // Running at Full Thrust, with Aft (yet)
    {
      Time     dt = a_ft     - IgnTime;
      fullMass    = FullMass - EngineMR * dt;
      fuelMass    = FuelMass - FuelMR   * dt;
      oxidMass    = OxidMass - OxidMR   * dt;
      absThrust   = ThrustVac;
      fuelMassDot = - FuelMR;
      oxidMassDot = - OxidMR;
    }
    else
    if (a_ft < CutOffTime)   // Running at Full Thrust, w/o Aft
    {
      Time dt0    = a_ft            - IgnTime;
      Time dt1    = a_ft            - AftJetTime;
      fullMass    = AftJetFullMass  - EngineMR * dt1;
      fuelMass    = FuelMass        - FuelMR   * dt0;
      oxidMass    = OxidMass        - OxidMR   * dt0;
      absThrust   = ThrustVac;
      fuelMassDot = - FuelMR;
      oxidMassDot = - OxidMR;
    }
    else                      // After the Engine Cut-Off: "Spent" Stage
    {
      fullMass    = CutOffFullMass;
      fuelMass    = CutOffFuelMass;
      oxidMass    = CutOffOxidMass;
      absThrust   = Force   (0.0);
      fuelMassDot = MassRate(0.0);
      oxidMassDot = MassRate(0.0);
    }

    // Verify the Masses:
    //  "egMass" includes EmptyMass and GasesMass:
    DEBUG_ONLY
    (
      Mass egMass   = fullMass - fuelMass - oxidMass;
      assert
        (egMass.ApproxEquals(a_ft < AftJetTime ? EGMassBefore : EGMassAfter));
    )
    // Also, the Fuel and Oxid masses must not be below the physical low
    // limits:
    assert(FuelRem <= fuelMass  && fuelMass <= FuelMass &&
           OxidRem <= oxidMass  && oxidMass <= OxidMass &&
           !(IsPos(fuelMassDot) || IsPos(oxidMassDot)));

    MassRate fullMassDot = fuelMassDot + oxidMassDot;

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
      Angle_deg  A   = a_chamber_defls[i];
      assert(Abs(A) <= GimbalAmpl);

      double sinA = Sin(To_Angle(A));
      double cosA = Cos(To_Angle(A));
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
    //-----------------------------------------------------------------------//
    // Moments of Inertia and the Center of Masses:                          //
    //-----------------------------------------------------------------------//
    // NB:
    // (*) "fuelLevel" and "oxidLevel" are hooks provided for testing purposes
    //     only;
    // (*) Gases (primarily He) are considered along with the EmptyMass; redist-
    //     ribution of Gases within Stage3 is NOT taken into account, whereas in
    //     reality, He is used for Fuel and Oxid tanks pressurisation, as fills
    //     the "empty" space above the remaining Fuel and Oxid.   However, this
    //     effect is considered negligible for Stage3 (though may have some imp-
    //     act on the Moments of Inertia for Stage2):
    //-----------------------------------------------------------------------//
    // Fuel: BulkME and Level:                                               //
    //-----------------------------------------------------------------------//
    assert(IsPos(fuelMass) && fuelMass <= FuelMass && FuelMass < FuelTankMC);
    Len fuelLevel(NaN<double>);

    ME fuelME =
      (fuelMass > FuelTankLowMidMC)
      ? // Fuel level is within the FuelTankUp:
        FuelTankUp .GetPropBulkME
          (tt, fuelMass - FuelTankLowMidMC,     fuelMassDot, &fuelLevel) +
        FuelLowMidME
      :
      (fuelMass > FuelTankLowMC)
      ? // Fuel level is within the FuelTankMid:
        FuelTankMid.GetPropBulkME
          (tt, fuelMass - FuelTankLowMC,        fuelMassDot, &fuelLevel) +
        FuelLowME
      :
        // Fuel level is within the FuelTankLow:
        FuelTankLow.GetPropBulkME(tt, fuelMass, fuelMassDot, &fuelLevel);
    assert(IsPos(fuelLevel));

    //-----------------------------------------------------------------------//
    // Oxid: BulkME and Level:                                               //
    //-----------------------------------------------------------------------//
    Len oxidLevel(NaN<double>);
    assert(IsPos(oxidMass) && oxidMass <= OxidMass && OxidMass < OxidTankMC);

    ME oxidME =
      (oxidMass > OxidTankLowMidMC)
      ? // Oxid level is within the OxidTankUp:
        OxidTankUp .GetPropBulkME
          (tt, oxidMass - OxidTankLowMidMC,     oxidMassDot, &oxidLevel) +
        OxidLowMidME
      :
      (oxidMass > OxidTankLowMC)
      ? // Oxid level is within the OxidTankMid:
        OxidTankMid.GetPropBulkME
          (tt, oxidMass - OxidTankLowMC,        oxidMassDot, &oxidLevel) +
        OxidLowME
      :
        // Oxid level is within the OxidTankLow:
        OxidTankLow.GetPropBulkME(tt, oxidMass, oxidMassDot, &oxidLevel);
    assert(IsPos(oxidLevel));

    //-----------------------------------------------------------------------//
    // Full = (Empty+Gases) + Fuel + Oxid:                                   //
    //-----------------------------------------------------------------------//
    ME fullME = (a_ft < AftJetTime ? EGBeforeME : EGAfterME) + fuelME + oxidME;

    // Double-check the Masses and the MassRate:
    assert(fullME.GetMass   ().ApproxEquals(fullMass)    &&
           fuelME.GetMass   ().ApproxEquals(fuelMass)    &&
           oxidME.GetMass   ().ApproxEquals(oxidMass)    &&
           fullME.GetMassDot().ApproxEquals(fullMassDot) &&
           fuelME.GetMassDot().ApproxEquals(fuelMassDot) &&
           oxidME.GetMassDot().ApproxEquals(oxidMassDot));

    //-----------------------------------------------------------------------//
    // Make the Result:                                                      //
    //-----------------------------------------------------------------------//
    return StageDynParams<LVSC::Soyuz21b>
    {
      .m_ft          = a_ft,
      .m_fullMass    = fullMass,
      .m_fuelMass    = fuelMass,
      .m_oxidMass    = oxidMass,
      .m_fullMassDot = fullMassDot,
      .m_com         = fullME.GetCoM    (),
      .m_comDots     = fullME.GetCoMDots(),
      .m_mois        = fullME.GetMoIs   (),
      .m_moiDots     = fullME.GetMoIDots(),
      // NB: BOTH ECOS TS and VecTS are equal to "tt" here:
      .m_thrust      = ForceVEmb<LVSC::Soyuz21b>
                         (tt, tt, thrustX, thrustY, thrustZ)
    };
  }
}
// End namespace SpaceBallistics
