// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/Soyuz-2.1b/Stages1-2.hpp":            //
//                      Mathematical Model of "Soyuz-2.1b"                   //
//     Stage2 (Block 'A') and Stage1 Boosters (Blocks 'B', 'V', 'G', 'D')    //
//===========================================================================//
#pragma  once
#ifdef   STAGE1
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"
#include "SpaceBallistics/PhysEffects/PlateAeroDyn.hpp"
#else
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#endif
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.h"
#include "SpaceBallistics/ME/TanksPressrn.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage2::GetDynParams":                                        //
  //=========================================================================//
# ifdef STAGE1
  template<char Block>
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage1_Booster<Block>::GetDynParams
# else
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage2::GetDynParams
# endif
  (
    FlightTime             a_ft,
    VernDeflections const& a_vern_defls,
    Pressure               a_p              // Curr Atmospheric Pressure
#   ifdef STAGE1
    ,
    // Extra params for AeroFin Ctls:
    Density                a_rho,           // Curr Atmospheric Density
    ME::VelVE const&       a_v,             // In the ECOS SnapShot
    ME::VelVE const&       a_w,             // ditto
    Angle_deg              a_aerofin_defl   // AeroFin deflection angle
#   endif
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    // We currently do not allow any times prior to LiftOff:
    assert(SC::LiftOffTime <= a_ft);

    // The absolute curr time must be defined in this case:
    TT  tt = TT(a_ft);
    assert(!tt.IsUnDef());

    // The atmospheric pressure is >= 0 obviously:
    assert(!IsNeg(a_p));

    // Check the MassRates at FullThrust:
    static_assert((OxidMR + FuelMR + H2O2MR).ApproxEquals(FullMR));

    // Check the Timings:
    static_assert(SC::LiftOffTime <= FullThrustTime &&
                  FullThrustTime  <  ThrottlTime    &&
                  ThrottlTime     <  MainCutOffTime &&
                  MainCutOffTime  <= CutOffTime);

    namespace EAM = EarthAtmosphereModel;

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//
    // XXX: We assume for simplicity that both "Vac" and "SL" Thrust values are
    // throttled in the same proportion, and the throttling levels  of the Main
    // Engines and the Verniers are the same:
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

    Time      dt0       = a_ft - SC::LiftOffTime;
    assert(!IsNeg(dt0));

#   ifdef STAGE1
    if (a_ft < FullThrustTime)
    {
      // The first seconds after LiftOff, for Stage1 only. For Stage2, this mode
      // is not present, since it develops FullThrust already at LiftOff time:
      fuelMass          = FuelMass0  - FuelMRL * dt0;
      oxidMass          = OxidMass0  - OxidMRL * dt0;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * dt0;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * dt0;

      fuelMassDot       = - FuelMRL;
      oxidMassDot       = - OxidMRL;
      h2o2MassDot       = - H2O2MR;    // H2O2 and LiqN2 rates are constant...
      liqN2MassDot      = - LiqN2MR;   //

      absMainThrustVac  = ThrustMainVac  * LiftOffThrottlLevel;
      absMainThrustSL   = ThrustMainSL   * LiftOffThrottlLevel;
      absVernThrustVac1 = ThrustVernVac1 * LiftOffThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * LiftOffThrottlLevel;
    }
    else
#   endif
    if (a_ft < ThrottlTime)
    {
      // Full-Thrust Mode, for both Stage1 and Stage2.
      // For Stage2, FullThrustTime == LiftOffTime:
      Time    dtF       = a_ft - FullThrustTime;

      fuelMass          = FuelMassF  - FuelMR  * dtF;
      oxidMass          = OxidMassF  - OxidMR  * dtF;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * dt0;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * dt0;

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
    if (a_ft < MainCutOffTime)
    {
      // Throttled Mode (for both Main Engine and Verniers), Stage1 and Stage2.
      // For Stage1, MainCutOffTime = CutOffTime  (ie the Main Engine  and the
      // Verniers are cut off simultaneously),   whereas for Stage2,
      // MainCutOffTime < CutOffTime:
      Time  dtT         = a_ft       - ThrottlTime;
      assert(!IsNeg(dtT));

      fuelMass          = FuelMassT  - FuelMRT * dtT;
      oxidMass          = OxidMassT  - OxidMRT * dtT;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * dt0;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * dt0;

      fuelMassDot       = - FuelMRT;
      oxidMassDot       = - OxidMRT;
      h2o2MassDot       = - H2O2MR;
      liqN2MassDot      = - LiqN2MR;

      absMainThrustVac  = ThrustMainVac  * ShutDownThrottlLevel;
      absMainThrustSL   = ThrustMainSL   * ShutDownThrottlLevel;
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
#   ifndef STAGE1
    else
    if (a_ft < CutOffTime)
    {
      // For Stage2 only: In this mode, only the Verniers continue operation,
      // in the Throttled Mode, prior to the complete CutOff:
      Time  dtM         = a_ft       - MainCutOffTime;
      assert(!IsNeg(dtM));

      fuelMass          = FuelMassM  - FuelMRM * dtM;
      oxidMass          = OxidMassM  - OxidMRM * dtM;
      h2o2Mass          = H2O2Mass0  - H2O2MR  * dt0;
      liqN2Mass         = LiqN2Mass0 - LiqN2MR * dt0;

      fuelMassDot       = - FuelMRM;
      oxidMassDot       = - OxidMRM;
      h2o2MassDot       = - H2O2MR;
      liqN2MassDot      = - LiqN2MR;

      absMainThrustVac  = Force   (0.0);
      absMainThrustSL   = Force   (0.0);
      absVernThrustVac1 = ThrustVernVac1 * ShutDownThrottlLevel;
      absVernThrustSL1  = ThrustVernSL1  * ShutDownThrottlLevel;
    }
#   endif
    else
    {
      // After CuttOff, for Stage1 and Stage2: Trivial:
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
    // H2O2 spent since LiftOff  (once again, N2O2 is NOT spent):
    Mass  fuelSpent     = FuelMass0 - fuelMass;
    Mass  oxidSpent     = OxidMass0 - oxidMass;
    Mass  h2o2Spent     = H2O2Mass0 - h2o2Mass;
    assert(!(IsNeg(fuelSpent) || IsNeg(oxidSpent) || IsNeg(h2o2Spent)));

    Mass  fullMass      = FullMass0 - (fuelSpent  + oxidSpent   + h2o2Spent);
    DEBUG_ONLY(MassRate fullMassDot = fuelMassDot + oxidMassDot + h2o2MassDot;)

    // Check the Curr Masses and MassRates:
    assert(MinEndMass <= fullMass  && fullMass  <= FullMass0  &&
           FuelRem    <= fuelMass  && fuelMass  <= FuelMass0  &&
           OxidRem    <= oxidMass  && oxidMass  <= OxidMass0  &&
           H2O2Rem    <= h2o2Mass  && h2o2Mass  <= H2O2Mass0  &&
           LiqN2Rem   <= liqN2Mass && liqN2Mass <= LiqN2Mass0 &&
           !(IsPos(fuelMassDot) || IsPos(oxidMassDot) || IsPos(h2o2MassDot) ||
             IsPos(liqN2MassDot)));

    // XXX: "h2o2Mass" is not returned to the CallER; it is only used in the
    // MoIs and CoM computations below...

    //-----------------------------------------------------------------------//
    // Thrust Vector:                                                        //
    //-----------------------------------------------------------------------//
    // Adjust the Thrust values to the curr atmospheric pressure:
    //
    double slC  = double(a_p / EAM::P0); // 1 when we are @ SL
    double vacC = 1.0 - slC;             // 1 when we are in Vac
    assert(0.0 <= slC && slC <= 1.0 && 0 <= vacC && vacC <= 1.0);

    Force absMainThrust  = vacC * absMainThrustVac  + slC * absMainThrustSL;
    Force absVernThrust1 = vacC * absVernThrustVac1 + slC * absVernThrustSL1;

    Force thrustX = absMainThrust;
    Force thrustY;  // Initially 0
    Force thrustZ;  // ditto

    // Consider VernGimbalAngles of all Vernier Chambers (2 for Stage1, 4 for
    // Stage2). The formulas here are the same as for Stage3, because the Stage1
    // and Stage2 Verniers and Stage3 Main (Gimbaled) Chambers are deflected in
    // the same planes!
    // The number of Verniers:
    int const NV = int(a_vern_defls.size());
#   ifdef STAGE1
    assert(NV == 2);
#   else
    assert(NV == 4);
#   endif

    for (int i = 0; i < NV; ++i)
    {
      Angle_deg  A   = a_vern_defls[size_t(i)];
      assert(Abs(A) <= VernGimbalAmpl);

      double sinA = Sin(To_Angle(A));
      double cosA = Cos(To_Angle(A));

      // The X component of the Vernier Thrust has the same expression for any
      // Venrnier and for both Stage1 and Stage2:
      thrustX    += absVernThrust1 * cosA;

      // The YZ component of the Vernier Thrust:
      // (*) the ThrustVector rotation in the YZ plane is OPPOSITE  to the
      //     corresp Vernier Deflection;
      // (*) for Stage1, the Vernier Deflection Plane is determined by the
      //     BlockID ('B', 'V', 'G', 'D'), and is the same for both Verniers
      //     (though their Deflection Angles may differ);
      // (*) for Stage2, the Vernier Deflection Plane is determined by the
      //     VernierID ("i"), and
      //     (i=0) <-> 'B', (i=1) <-> 'V', (i=2) <-> 'G', (i=3) <-> 'D':
#     ifdef STAGE1
      switch (Block)
#     else
      switch (i)
#     endif
      {
      case  0 :
      case 'B':
        thrustZ += absVernThrust1 * sinA;
        break;
      case  1 :
      case 'V':
        thrustY -= absVernThrust1 * sinA;
        break;
      case  2 :
      case 'G':
        thrustZ -= absVernThrust1 * sinA;
        break;
      case  3 :
      case 'D':
        thrustY += absVernThrust1 * sinA;
        break;
      default:
        assert(false);
      }
    }
    //-----------------------------------------------------------------------//
    // Moments of Inertia and the Center of Masses:                          //
    //-----------------------------------------------------------------------//
    // XXX: Similar to Stage3, "fuelLevel", "oxidLevel" and "h2o2Level" are
    // provided for testing purposes only:
    //-----------------------------------------------------------------------//
    // Fuel: BulkME and Level, for both Stage1 and Stage2:                   //
    //-----------------------------------------------------------------------//
    Len fuelLevel(NaN<double>);

    ME fuelME =
      (fuelMass > FuelTankBtmMidMC)
      ? // Fuel level is within the FuelTankTop:
        FuelTankTop.GetPropBulkME
          (tt, fuelMass - FuelTankBtmMidMC,     fuelMassDot,  &fuelLevel) +
        FuelBtmMidME
      :
      (fuelMass > FuelTankBtmMC)
      ? // Fuel level is within the FuelTankMid:
        FuelTankMid.GetPropBulkME
          (tt, fuelMass - FuelTankBtmMC,        fuelMassDot,  &fuelLevel) +
        FuelBtmME
      :
        // Fuel level is within the FuelTankBtm:
        FuelTankBtm.GetPropBulkME(tt, fuelMass, fuelMassDot,  &fuelLevel);

    assert(IsPos(fuelLevel));

    //-----------------------------------------------------------------------//
    // Oxid: BulkME and Level, for both Stage1 and Stage2:                   //
    //-----------------------------------------------------------------------//
    Len oxidLevel(NaN<double>);

    ME oxidME =
      (oxidMass > OxidTankBtmLowUpMC)
      ? // Oxid level is within the OxidTankTop:
        OxidTankTop.GetPropBulkME
          (tt, oxidMass - OxidTankBtmLowUpMC,   oxidMassDot,  &oxidLevel) +
        OxidBtmLowUpME
      :
      (oxidMass > OxidTankBtmLowMC)
      ? // Oxid level is within the OxidTankUp:
        OxidTankUp .GetPropBulkME
          (tt, oxidMass - OxidTankBtmLowMC,     oxidMassDot,  &oxidLevel) +
        OxidBtmLowME
      :
      (oxidMass > OxidTankBtmMC)
      ? // Oxid level is within the OxidTankLow:
        OxidTankLow.GetPropBulkME
          (tt, oxidMass - OxidTankBtmMC,        oxidMassDot,  &oxidLevel) +
        OxidBtmME
      :
        // Oxid level is within the OxidTankBtm:
        OxidTankBtm.GetPropBulkME(tt, oxidMass, oxidMassDot,  &oxidLevel);

    assert(IsPos(oxidLevel));

    //-----------------------------------------------------------------------//
    // H2O2: BulkME and Level:                                               //
    //-----------------------------------------------------------------------//
    // Stage2 has an extra H2O2 Tank section ("Mid") which is absent in Stage1:
    //
    Len h2o2Level(NaN<double>);

    ME h2o2ME =
#     ifndef STAGE1
      (h2o2Mass > H2O2TankBtmMidMC)
      ? // H2O2 level is within the H2O2TankTop (Stage1 only):
        H2O2TankTop.GetPropBulkME
          (tt, h2o2Mass - H2O2TankBtmMidMC,     h2o2MassDot,  &h2o2Level) +
        H2O2BtmMidME
      :
#     endif
      (h2o2Mass > H2O2TankBtmMC)
      ? // H2O2 level is within the H2H2TankMid(Stage2) or H2H2TankTop(Stage1):
#     ifndef STAGE1
        H2O2TankMid.GetPropBulkME
#     else
        H2O2TankTop.GetPropBulkME
#     endif
          (tt, h2o2Mass - H2O2TankBtmMC,        h2o2MassDot,  &h2o2Level) +
        H2O2BtmME
      :
        // H2O2 level is within the H2O2TankBtm:
        H2O2TankBtm.GetPropBulkME(tt, h2o2Mass, h2o2MassDot,  &h2o2Level);

    assert(IsPos(h2o2Level));

    //-----------------------------------------------------------------------//
    // LiqN2: BulkME and Level, for both Stage1 and Stage2:                  //
    //-----------------------------------------------------------------------//
    Len liqN2Level(NaN<double>);

    ME liqN2ME =
      (liqN2Mass > LiqN2TankBtmMC)
      ? // LiqN2 level is within the LiqN2TankTop:
        LiqN2TankTop.GetPropBulkME
          (tt, liqN2Mass - LiqN2TankBtmMC,        liqN2MassDot, &liqN2Level) +
        LiqN2BtmME
      :
        // LiqN2 level is within the LiqN2TankBtm:
        LiqN2TankBtm.GetPropBulkME(tt, liqN2Mass, liqN2MassDot, &liqN2Level);
    assert(IsPos(liqN2Level));

    //-----------------------------------------------------------------------//
    // Gaseous N2:                                                           //
    //-----------------------------------------------------------------------//
    // Unlike Stage3, here we cannot neglect the effects  of the Tank Pressuri-
    // sation Gases (in this case, N2), because their mass is relatively large.
    // N2 is NOT spent (exhausted), but its distribution within the Tanks chan-
    // ges over the time. Gasesous N2 is concentrated in the upper parts of the
    // Fuel and Oxid Tanks, gradually increasing in volume.
    // The following logic is the same for Stage1 and Stage2:
    //
    // The "complementary" ("CME") parts of the Fuel, Oxid and LiqN2 bulks
    // (which are actually filled with Gaseous N2):
    //
    Mass     gasN2Mass  = N2Mass - liqN2Mass;
    assert  (GasN2Mass <= gasN2Mass && gasN2Mass < N2Mass - LiqN2Rem);

    MassRate gasN2MassDot = - liqN2ME.GetMassDot();
    assert  (!IsNeg(gasN2MassDot));

    ME gasN2ME =
      TanksPressrn<LVSC::Soyuz21b>
      (
        FuelME,    fuelME,
        OxidME,    oxidME,
        gasN2Mass, gasN2MassDot
      );
    assert((liqN2ME.GetMass() + gasN2ME.GetMass()).ApproxEquals(N2Mass));

    //-----------------------------------------------------------------------//
    // We can now construct the FullME (for both Stage1 and Stage2):         //
    //-----------------------------------------------------------------------//
    ME fullME = EmptyME + fuelME + oxidME + h2o2ME + liqN2ME + gasN2ME;

    // Double-check the Masses and the MassRate:
    assert(fuelME .GetMass   ().ApproxEquals(fuelMass));
    assert(oxidME .GetMass   ().ApproxEquals(oxidMass));
    assert(h2o2ME .GetMass   ().ApproxEquals(h2o2Mass));
    assert(liqN2ME.GetMass   ().ApproxEquals(liqN2Mass));
    assert(gasN2ME.GetMass   ().ApproxEquals(gasN2Mass));
    assert(fullME .GetMass   ().ApproxEquals(fullMass));
    assert(fuelME .GetMassDot().ApproxEquals(fuelMassDot));
    assert(oxidME .GetMassDot().ApproxEquals(oxidMassDot));
    assert(h2o2ME .GetMassDot().ApproxEquals(h2o2MassDot));
    assert(liqN2ME.GetMassDot().ApproxEquals(liqN2MassDot));
    assert(fullME .GetMassDot().ApproxEquals(fullMassDot));

    //-----------------------------------------------------------------------//
    // Moments from AeroFins:                                                //
    //-----------------------------------------------------------------------//
#   ifdef STAGE1
    // Normal Unit Vector to the AeroFin:
    DimLessVEmb<LVSC::Soyuz21b> n;   // Initially, all components are 0s

    Angle  theta = To_Angle(a_aerofin_defl);
    n.x()        = Cos(theta);

    switch (Block)
    {
      case 'B':
        n.z() =   Sin(theta);
        break;
      case 'V':
        n.y() = - Sin(theta);
        break;
      case 'G':
        n.z() = - Sin(theta);
        break;
      case 'D':
        n.y() =   Sin(theta);
        break;
      default:
        assert(false);
    }
    // AeroDynamic Force  on the AeroFin:
    ForceVEmb<LVSC::Soyuz21b> aeroFinF =
      PlateAeroDyn(a_v, a_w, n, a_p, a_rho, EAM::GammaAir, AeroFinS);
#   endif

    //-----------------------------------------------------------------------//
    // Make the Result:                                                      //
    //-----------------------------------------------------------------------//
    return StageDynParams<LVSC::Soyuz21b>
    {
      .m_ft          = a_ft,
      .m_fullMass    = fullMass,
      .m_fuelMass    = fuelMass,
      .m_oxidMass    = oxidMass,
      .m_fullMassDot = fullME.GetMassDot(),
      .m_com         = fullME.GetCoM    (),
      .m_comDots     = fullME.GetCoMDots(),
      .m_mois        = fullME.GetMoIs   (),
      .m_moiDots     = fullME.GetMoIDots(),
      .m_thrust      = ForceVEmb<LVSC::Soyuz21b>(tt, thrustX, thrustY, thrustZ)
    };
  }
}
// End namespace SpaceBallistics
