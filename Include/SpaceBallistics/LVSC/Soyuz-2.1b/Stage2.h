// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h":                //
//         Mathematical Model of the "Soyuz-2.1b" Stage2 ("Block A")         //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/ME/TrConeSpherSegm.hpp"
#include "SpaceBallistics/ME/ToricSegms.hpp"
#include "SpaceBallistics/LVSC/Propellants.h"
#include "SpaceBallistics/LVSC/StageDynParams.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Consts.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Head.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC  = Soyuz21b_Consts;
  namespace SH  = Soyuz21b_Head;

  // All "MechElements" are instantiated with "LVSC::Soyuz21b":
  using     ME   = MechElement   <LVSC::Soyuz21b>;
  using     PM   = PointMass     <LVSC::Soyuz21b>;
  using     TrC  = TrCone        <LVSC::Soyuz21b>;
  using     SpS  = SpherSegm     <LVSC::Soyuz21b>;
  using     Tor  = ToricSegm     <LVSC::Soyuz21b>;
  using     DCyl = DoubleCylinder<LVSC::Soyuz21b>;
  using     S3   = Soyuz21b_Stage3;

  //=========================================================================//
  // "Soyuz21b_Stage2" Class:                                                //
  //=========================================================================//
  // Similar to "Soyuz21b_Stage3", this class is actually a namespace, with most
  // members being "constexpr" "static"s:
  //
  class Soyuz21b_Stage2
  {
  private:
    //=======================================================================//
    // No objects can be created of this class:                              //
    //=======================================================================//
    Soyuz21b_Stage2() = delete;

  public:
    //=======================================================================//
    // Geometry:                                                             //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // InterStageGrid, DeflectorCone, EquipmentBay:                          //
    //-----------------------------------------------------------------------//
    // Stage2 Top is the Bottom of Stage3. It is also the top of the InterStage
    // Grid:
    constexpr static Len      TopX          = SC::X0 - SC::Stage3Len;

    // The IntegrStageGrid:
    constexpr static Len      GridH         = 0.975_m;

    // DeflectorCone Height and Base Diameter. The Height is < GridH, so the
    // DeflectorCone is within the Grid:
    constexpr static Len      DeflConeH     = GridH / 2.0;
    constexpr static Len      TopD          = 2.700_m;  // Slightly > Stage3D

    // The Grid musy be tall enough to accommodate the Stage3 Nozzles extending
    // from the Stage3 itself, and the DeflCone:
    static_assert(GridH > S3::NozzlesExtH + DeflConeH);

    // Equipment Bay: Cylindrical or expanding TrCone? Information differs in
    // different drawings. Assume Cylindrical, of the same "TopD" diameter:
    constexpr static Len      EquipBayH     = 0.595_m;

    //-----------------------------------------------------------------------//
    // OxidTank:                                                             //
    //-----------------------------------------------------------------------//
    // The over-all Upper (Expanding) TrCone: Contains the Top SpherSegm of the
    // OxidTank, and the Upper TrCone of the OxidTank itself. The diameter inc-
    // reases from "TopD" to "MaxD":
    constexpr static Len      UpperTrCH     = 5.840_m;  // Or 5.835 m?
    constexpr static Len      MaxD          = 2.950_m;  // Max for entire Stage2

    // The over-all Middle (Narrowing-Down) TrCone: Contains the Low TrCone  of
    // the OxidTank. The diameter decreses from MaxD to MinD. The Middle TrCone
    // is also the top junction point of the Stage1 StrapOn Boosters, wo we wil
    // need its angle:
    constexpr static Len      MidTrCH       = 4.55_m;
    constexpr static Len      MinD          = 2.05_m;
    constexpr static double   TanAlphaMid   =
      double((MaxD - MinD) / 2.0 / MidTrCH);

    // Top OxidTank SpherSegm. Its top-most point (the pole) is at the level of
    // the notional bottom base of the EquipBay:
    constexpr static Len      OxidTankTopH  = 0.8_m;    // Approx
    constexpr static Len      OxidTankTopD  =
      TopD + double(OxidTankTopH / UpperTrCH) * (MaxD - TopD);

    // Height of the Upper TrCone part of the OxidTank: extends to the end of
    // the over-all  Upper TrCone:
    constexpr static Len      OxidTankUpH   = UpperTrCH - OxidTankTopH;

    // Height of the Lower TrCone part of the OxidTank (extends ALMOST to the
    // end of the Middle TrCone, so its diameter is slightly > MinD):
    constexpr static Len      OxidTankBtmH  = 0.5_m;    // Approx
    constexpr static Len      OxidTankBtmH2 = OxidTankBtmH / 2.0;
    constexpr static Len      OxidTankLowH  = MidTrCH   - OxidTankBtmH2;
    constexpr static Len      OxidTankLowD  =
      MinD + (MaxD - MinD) * double(OxidTankBtmH2 / MidTrCH);

    // Bottom SpherSegm of OxidTank: D=OxidTankLowD, H=OxidTankBtmH

    //-----------------------------------------------------------------------//
    // FuelTank:                                                             //
    //-----------------------------------------------------------------------//
    // Located entirely in the Cylindrical Low section of Stage2. The Diameter
    // of this section is "MinD", and over-all length is:
    constexpr static Len      LowCylH       = 15.178_m;

    // There is a gap between the OxidTankBtm and FuelTankTop. Its exact size
    // is not known, so assume something similar to Stage3:
    constexpr static Len      OxidFuelGap   = 0.15_m;

    // Top SpherSegm:
    constexpr static Len      FuelTankTopH  = OxidTankBtmH;   // XXX: Approx...

    // Bottom SpherSegm:
    constexpr static Len      FuelTankBtmH  = FuelTankTopH;

    // TrCone Enclosure of the top-half of "OxidTankBtm" (NB: the lower-half is
    // enclosed by the "OxidFuelEnc", see below):
    constexpr static Len      OxidBtmEnclH  = OxidTankBtmH2;

    // Length of the main cylindrical section of the FuelTank: Known only
    // approximately. We assume (see above) that half of the OxidTankBtm secti-
    // on is in the "MidTrCone", and the "FuelTankBtm" SpherSegm is beyond the
    // 10.378_m quoted below:
    constexpr static Len      OxidFuelEnclH =
      OxidTankBtmH2 + OxidFuelGap + FuelTankTopH;

    constexpr static Len      FuelTankMidH  = 10.378_m - OxidFuelEnclH;

    //-----------------------------------------------------------------------//
    // H2O2 and N2 Tanks:                                                    //
    //-----------------------------------------------------------------------//
    // Top and Bottom parts are Semi-ToricSegms, the Mid part is a Double
    // Cylinder. All sizes are approximate:
    constexpr static Len      H2O2TankOutR   = MinD / 2.0;
    constexpr static Len      H2O2TankInR    = 0.63_m;
    constexpr static Len      H2O2TankMinoR  =
      (H2O2TankOutR - H2O2TankInR) / 2.0;
    constexpr static Len      H2O2TankCylH   = 1.2_m;

    // N2 Tank is a Full (Circular) Torus:
    constexpr static Len      LiqN2TankMinoR = 0.2_m;

    // We assume that the two Gaps (FuelTank-H2O2Tank and H2O2Tank-LiqN2Tank)
    // are the same, and derive them from the over-all compartment length. In
    // that part of the Tail compartment, the following are located:
    //   FuelTankBtm, H2O2Gap, H2O2Tank, H2O2Gap, LiqN2Tank:
    // Make this gap pretty small, otherwise the Engine may not fit:
    constexpr static Len      H2O2Gap        = 0.05_m;

    //-----------------------------------------------------------------------//
    // Over-All length of Stage2:                                            //
    //-----------------------------------------------------------------------//
    // (From the top of the InterStageGrid to the end of the LowCyl, but w/o
    // the extending Engine Nozzles):
    constexpr static Len      H =
      GridH + EquipBayH + UpperTrCH + MidTrCH + LowCylH;
    static_assert(H.ApproxEquals(27.138_m));

    // Nozzles extension beyond "H".  XXX: This val is for the Verniers; the
    // Main Engine Nozzles extend slightly less, but for the moment, we disre-
    // gard the difference:
    constexpr static Len      NozzlesExtH   = 0.66_m;

    // RD-108A Height:
    constexpr static Len      EngineH       = 2.865_m;

    // The Lowest Co-Ord (the Nozzles End) which will be used to compute the
    // Thrust Momenta:
    constexpr static Len      NozzlesLowX   = TopX - H - NozzlesExtH;

    //=======================================================================//
    // Masses:                                                               //
    //=======================================================================//
    // EmptyMass: XXX: StarSem says 6545 kg:
    constexpr static Mass     EmptyMass     = 6450.0_kg;

    // Mass of the RD-108A engine (part of EmptyMass):
    constexpr static Mass     EngMass       = 1075.0_kg;

    // Masses of Fuel (Naftil) and Oxidiser. As for Stage3, FuelMass includes
    // extra 0.2% for the antifreeze (2-EtoxyEthanol):
    constexpr static Mass     FuelMass      = 26794.0_kg * 1.002;
                                                          // StarSem: 26300 (T1)
    constexpr static Mass     OxidMass      = 63709.0_kg; // StarSem: 63800
    constexpr static Mass     H2O2Mass      = 2636.0_kg;
    constexpr static Mass     N2Mass        = 509.0_kg;
    constexpr static Mass     FullMass      = EmptyMass + FuelMass + OxidMass +
                                              H2O2Mass  + N2Mass;
    // N2Mass is split into a Liquid and Gaseous Phases, the initial masses are:
    constexpr static Mass     LiqN2Mass     = 485.0_kg;
    constexpr static Mass     GasN2Mass     =  24.0_kg;
    static_assert(LiqN2Mass + GasN2Mass == N2Mass);

    // UnSpendable Remnants of the Fuel and Oxidiser in Stage2 at the engine
    // cut-off time.   NB: The Remnants are about 1% of the corresp initial
    // masses. NB: They are Technically UnSpendable, so they do NOT include the
    // "guarantee margins"(???):
    constexpr static Mass     FuelRem       = 272.0_kg * 1.002;
    constexpr static Mass     OxidRem       = 678.0_kg;
    constexpr static Mass     H2O2Rem       = 263.0_kg;
    constexpr static Mass     LiqN2Rem      = 60.0_kg;

    //=======================================================================//
    // RD-108A (14D21) Engine Performance:                                   //
    //=======================================================================//
    // Isp (SL/Vac, sec), assuming this is for the Engine as a whole:
    // 257.7/320.6 (EnergoMash, RosCosmos.RU/2115),
    // 257.5/320.4 (LPRE.DE),
    // 255.0/319.0 (StarSem;    this looks slightly low),
    // 268.1/333.6 (Kravchenko; this looks a way too high!);
    // Vac/SL ration is ~ 1.244 (except Kravchenko which is a clear outlier);
    //
    // Thrust (SL/Vac, tf):
    // 70/87       (RosCosmos.RU/2115: Main Chambers only, w/o the Verniers!)
    // 80.8/94     (Kravchenko with a typo corrected, EnergoMash, Wikipedia;
    //              the Vac value seems to grossly incorrect!)
    // 80.8/101.0  (LPRE.DE, StarSem, Zhelesnyakov).
    //
    // For both SL and Vac conditions, the Engine MassRate = Thrust/Isp must be
    // constant. Based on the above data,  the MassRate is indeed nearly const-
    // ant; it varies from 313.5 (SL) to 315.0 (Vac) kg/sec. We assume the lower
    // value (to extend the MaxBurnTime) which corresponds to SL:
    constexpr static Time     IspEngSL       = 257.7_sec;
    constexpr static Time     IspEngVac      = 320.6_sec;
    constexpr static Force    ThrustEngSL    = 80800.0_kg  * g0;
    constexpr static MassRate EngineMR       = ThrustEngSL / (IspEngSL * g0);
    // EngineMR = 313.54 kg/sec

    // Then the Vac vals: we can exactly match either Thrust or Isp, but not
    // both; so match the Isp:
    constexpr static Force    ThrustEngVac   = IspEngVac   * g0 * EngineMR;
    // ThrustEngVac = 100.52 tf, sufficiently close to 101 tf quoted above!

    // Now the separate Isp and Thrust vals for the Main Engine  and the Verni-
    // ers.   Starting from the MainThrust  of 70/87 and the above ThrustEng of
    // 80.8/100.52,  we get the VernThrust4 of 10.8/13.52.
    // IMPORTANT:
    // Vernier Chambers of Stage2 are installed in positions at +-Y, +-Z  axes,
    // ie towards the corresp Stage1 blocks. They are rotatable to +-45 degress
    // in the planes parallel to the X axis, ie, the Vernier Chambers installed
    // at +-Y are rotatable in the XZ plane, and those installed at +-Z are ro-
    // tatable in the XY plane:
    // We accept the following Isp vals for the Verniers (LPRE.DE):
    //
    constexpr static Time     IspVernSL      = 251.9_sec;
    constexpr static Time     IspVernVac     = 313.1_sec;

    // Then the VernMR4 (SL/Vac) is 42.9/43.2 kg/sec, which is reasonably close
    // to a constant. Thus, we select
    constexpr static MassRate VernMR4        = MassRate(43.0);
    constexpr static MassRate VernMR1        = VernMR4 / 4.0;
    // 10.75 kg/sec, short of 4.15+8.55=12.7 kg/sec given by Kravchenko...
    // Thus,
    constexpr static Force    ThrustVernSL4  = VernMR4 * g0 * IspVernSL;
    constexpr static Force    ThrustVernVac4 = VernMR4 * g0 * IspVernVac;
    // 10.83/13.46 tf, very close to the above estimates!

    constexpr static Force    ThrustVernSL1  = ThrustVernSL4  / 4.0;
    constexpr static Force    ThrustVernVac1 = ThrustVernVac4 / 4.0;
    // 2.708/3.366 tf,  still short of 35 kN = 3.57 tf given by StarSem...

    // We can now derive the MassRates and the Isp vals for the Main Engine:
    constexpr static MassRate MainMR         = EngineMR     - VernMR4;
    constexpr static Force    ThrustMainSL   = ThrustEngSL  - ThrustVernSL4;
    constexpr static Force    ThrustMainVac  = ThrustEngVac - ThrustVernVac4;
    // For the ThrustMain, we get 69.97/87.06 tf, very close to 70/87 tf quoted
    // by RosCosmos!

    // Separate Fuel and Oxid Rates for the whole Engine, obtained using the
    // Oxid/Fuel Ratio which we derive from the over-all Fuel and Oxid spend-
    // able masses:
    constexpr static double OxidFuelRat =
      double((OxidMass - OxidRem) / (FuelMass - FuelRem));
    // ~2.38, very close to 2.39 often quoted for RD-108A

    constexpr static double FuelPart = 1.0  / (1.0 + OxidFuelRat);
    constexpr static double OxidPart = OxidFuelRat * FuelPart;

    constexpr static MassRate FuelMR = EngineMR * FuelPart;
    constexpr static MassRate OxidMR = EngineMR * OxidPart;

    //-----------------------------------------------------------------------//
    // RD-108A Ignition Sequence:                                            //
    //-----------------------------------------------------------------------//
    using FT = FlightTime;

    // For RD-108A, we assume that FullThrust instant is the same as LiftOff
    // (Contact Separation), ie t0=0:
    constexpr static FT       FullThrustTime  = SC::LiftOffTime;

    // Preliminary thrust level (assumed to be a 25%) occurs notionally at -15
    // sec, for both the Main Engine and the Vernier Engines (XXX: in reality,
    // the ignition sequence is more complex).
    //
    // So the Masses at LiftOff (t0=0) are following. NB: H2O2 is NOT spent at
    // that point yet, the TurboPumps become operational only at t0=0; the init-
    // ial Propellants flow is autonomous:
    //
    constexpr static Time     IgnAdvance      = 15.0_sec;
    constexpr static double   IgnThrottlLevel = 0.25;
    constexpr static Mass     FuelMass0       =
      FuelMass - FuelMR * IgnAdvance * IgnThrottlLevel;
    constexpr static Mass     OxidMass0       =
      OxidMass - OxidMR * IgnAdvance * IgnThrottlLevel;

    // H2O2 and LiqN2 masses are unchanged prior to LiftOff:
    constexpr static Mass     H2O2Mass0       = H2O2Mass;
    constexpr static Mass     LiqN2Mass0      = LiqN2Mass;
    constexpr static Mass     GasN2Mass0      = GasN2Mass;

    // Thus, the FullMass at LiftOff (t0=0):
    constexpr static Mass     FullMass0       =
      FullMass - ((FuelMass + OxidMass) - (FuelMass0 + OxidMass0));
    static_assert(FullMass0 < FullMass);

    // This is the same as masses at FullThrustTime (tF=t0=0):
    constexpr static Mass     FuelMassF       = FuelMass0;
    constexpr static Mass     OxidMassF       = OxidMass0;

    //-----------------------------------------------------------------------//
    // RD-108A Shut-Down Sequence:                                           //
    //-----------------------------------------------------------------------//
    // XXX: The exact sequence for RD-108A ShutDown is unknown, and almost cer-
    // tainly differs from that of RD-108,  because the latter  would  imply a
    // much lower BurnTime. The following fugures are manually  adjusted to be
    // consistent with "CutOffTime" and "EngineMR". For simplicity  we  assume
    // that the Main and Vernier Chambers are throttled at the same time,   to
    // the same 50% level:
    constexpr static double   ShutDownThrottlLevel  = 0.50;
    constexpr static Time     ThrottlAdvance        = 6.0_sec;
    constexpr static Time     MainSDAdvance         = 5.0_sec;

    constexpr static FT       CutOffTime      = SC::Stage2CutOffTime;
    constexpr static FT       ThrottlTime     = CutOffTime - ThrottlAdvance;
    constexpr static FT       MainCutOffTime  = CutOffTime - MainSDAdvance;

    // Thus, the Propellant Mass spent during the ShutDown sequence is:
    constexpr static Mass     ShutDownSpentPropMass =
       ShutDownThrottlLevel *
        ((ThrottlAdvance - MainSDAdvance) * EngineMR +
         MainSDAdvance                    * VernMR4);

    // Thus, the Propellant Mass left for the FullThrust Mode:
    constexpr static Mass     FullThrustPropMass    =
      FuelMass0 + OxidMass0 - ShutDownSpentPropMass - FuelRem - OxidRem;
    static_assert(IsPos(FullThrustPropMass));

    // Fuel and Oxid Masses @ "ThrottlTime":
    constexpr static Mass     FuelMassT =
      FuelMass0 - FuelMR  * (ThrottlTime - SC::LiftOffTime);
    constexpr static Mass     OxidMassT =
      OxidMass0 - OxidMR  * (ThrottlTime - SC::LiftOffTime);

    // Between "ThrottlTime" and "MainCutOffTime", the following MassRates
    // apply:
    constexpr static MassRate FuelMRT   = FuelMR * ShutDownThrottlLevel;
    constexpr static MassRate OxidMRT   = OxidMR * ShutDownThrottlLevel;

    // Propellant Masses at the "MainCutOffTime":
    constexpr static Mass     FuelMassM =
      FuelMassT - FuelMRT * (MainCutOffTime - ThrottlTime);
    constexpr static Mass     OxidMassM =
      OxidMassT - OxidMRT * (MainCutOffTime - ThrottlTime);

    // Between "MainCutOffTime" and "CutOffTime", the following MassRates
    // apply (only the Verniers continue to burn, in the Throttled mode):
    constexpr static MassRate FuelMRM   =
      VernMR4 * ShutDownThrottlLevel * FuelPart;
    constexpr static MassRate OxidMRM   =
      VernMR4 * ShutDownThrottlLevel * OxidPart;

    //-----------------------------------------------------------------------//
    // RD-108A BurnTimes and Rates:                                          //
    //-----------------------------------------------------------------------//
    // Then the Max Duration at FullThrust is:
    constexpr static Time  MaxFullThrustDur  =
      FullThrustPropMass / EngineMR;

    // Then the Max Time of RD-108A operation (to CutOff) is the following:
    constexpr static Time  MaxBurnDur  = MaxFullThrustDur + ThrottlAdvance;
    constexpr static FT    MaxBurnTime = SC::LiftOffTime  + MaxBurnDur;

    // "CutOffTime" must be less than (but close to) the "MaxBurnTime":
    static_assert(CutOffTime < MaxBurnTime);

    // NB: In addition, there is a MassRate due to H2O2 burning   to drive the
    // TurboPumps; however, this MassRate does not formally participate in the
    // Thrust. We assume that H2O2 flow begins at FullThrust / LiftOff time
    // and lasts for the whole "MaxBurnDur", with an exactly constant rate:
    //
    constexpr static MassRate H2O2MR  = (H2O2Mass - H2O2Rem) / MaxBurnDur;
    static_assert(IsPos(H2O2MR));

    // So  the total MassRate (at Full Thrust) is:
    constexpr static MassRate FullMR  = EngineMR + H2O2MR;

    // We assume that vaporisation of LiqN2 is also exactly  linear over time,
    // from FullThrust/LiftOff to the MaxBurnTime, and that initially, all of
    // N2 is in the Liquid form (XXX: which is not exactly correct).  However,
    // unlike H2O2, N2 is not exhaused, only re-distributed over the Tank vol-
    // umes becoming available:
    //
    constexpr static MassRate LiqN2MR = (LiqN2Mass - LiqN2Rem) / MaxBurnDur;
    static_assert(IsPos(LiqN2MR));

    // For testing: Minimal Mass of the Spent Stage2 (with all Remnants at their
    // minimal physical levels):
    constexpr static Mass     MinEndMass =
      EmptyMass + FuelRem + OxidRem + H2O2Rem + N2Mass;

    //-----------------------------------------------------------------------//
    // Fuel, Oxid, H2O2 and LiqN2 Masses @ "CutOffTime":                     //
    //-----------------------------------------------------------------------//
    constexpr static Mass    FuelMassC  =
      FuelMassM - FuelMRM * (CutOffTime - MainCutOffTime);
    constexpr static Mass    OxidMassC  =
      OxidMassM - OxidMRM * (CutOffTime - MainCutOffTime);
    constexpr static Mass    H2O2MassC  =
      H2O2Mass  - H2O2MR  * (CutOffTime - SC::LiftOffTime);
    constexpr static Mass    LiqN2MassC =
      LiqN2Mass - LiqN2MR * (CutOffTime - SC::LiftOffTime);

    // Checks:
    static_assert
      (FuelMassC  > FuelRem && OxidMassC > OxidRem && H2O2MassC > H2O2Rem &&
       LiqN2MassC > LiqN2Rem);

    static_assert((FuelMassC + OxidMassC).ApproxEquals
                  (FuelMassT + OxidMassT - ShutDownSpentPropMass));

  private:
    //=======================================================================//
    // "MechElement"s: "Proto"s with Yet-UnKnown Masses:                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Top-Level ME Protos:                                                  //
    //-----------------------------------------------------------------------//
    // XXX: Inter-Stage Grid (between Stages 2 and 3)  is  not modeled at all.
    // Its mass is included into the EmptyMass of Stage2, and will be "spread"
    // over the masses of all MEs defined below.  The reason is that the Grid
    // cannot (as yet) be represented as any supported "ME".  This will have a
    // minor effect on MoIs accuracy...
    //
    // "DeflectorCone" at the top of Stage2. The Ctor requires the vertex point
    // (not the base!):
    constexpr static TrC DeflectorConeProto =
      TrC
      (
        TT::UnDef(),
        TopX - GridH + DeflConeH,       // Cone vertex point
        0.0_m,                          //  ... so its D=0
        TopD,                           // Lower Base  D
        DeflConeH,
        Density(0.0)                    // No Propellant there
      );

    // "EquipBay" is assumed to be a Cylinder, its UpperBase is the LowBase of
    // "DeflectorCone":
    constexpr static TrC EquipBayProto =
      TrC
      (
        TT::UnDef(),
        DeflectorConeProto.GetLow()[0],
        TopD,
        EquipBayH,
        Density(0.0)                    // No propellant there
      );
    static_assert(EquipBayProto.GetUp()[0].ApproxEquals(TopX - GridH));

    //-----------------------------------------------------------------------//
    // OxidTank Protos:                                                      //
    //-----------------------------------------------------------------------//
    // The Enclosure of the OxidTankTop. Its UpperBase is the LowBase of the
    // "EquipBay", so its D is "TopD" as well:
    //
    constexpr static TrC  OxidTankTopEnclProto =
      TrC
      (
        TT::UnDef(),
        EquipBayProto.GetLow()[0],
        TopD,
        OxidTankTopD,
        OxidTankTopH,
        Density(0.0)                    // No propellant in it by itself...
      );

    // The Top SpherSegm of the OxidTank. Its UpperBase is the LowBase of the
    // above "OxidTankTopEnclProto":
    //
    constexpr static SpS OxidTankTopProto =
      SpS
      (
        TT::UnDef(),
        true,                           // Facing Up
        OxidTankTopEnclProto.GetLow()[0],
        OxidTankTopD,
        OxidTankTopH,
        Propellants::LOxDens
      );

    // The Upper TrCone of the OxidTank. Its UpperBase is the LowBase of the
    // "OxidTankTopProto" (same as the LowBase of "OxidTankTopEnclProto"):
    //
    static_assert(OxidTankTopEnclProto.GetLow()[0] ==
                  OxidTankTopProto    .GetLow()[0]);

    constexpr static TrC OxidTankUpProto =
      TrC
      (
        TT::UnDef(),
        OxidTankTopProto.GetLow()[0],
        OxidTankTopD,
        MaxD,
        OxidTankUpH,
        Propellants::LOxDens
      );

    // The Lower TrCone of the OxidTank. Stretches ALMOST to the end of the
    // "MidTrCone":
    constexpr static TrC OxidTankLowProto =
      TrC
      (
        TT::UnDef(),
        OxidTankUpProto.GetLow()[0],
        MaxD,
        OxidTankLowD,
        OxidTankLowH,
        Propellants::LOxDens
      );

    // TrCone Enclosure of the top-half of the "OxidTankBtm":
    constexpr static TrC OxidTankBtmEnclProto =
      TrC
      (
        TT::UnDef(),
        OxidTankLowProto.GetLow()[0],
        OxidTankLowD,
        MinD,
        OxidTankBtmH2,
        Density(0.0)                // No propellant there
      );

    // The Bottom SpherSegm of the OxidTank. Its Base is the LowBase of the
    // "OxidTankLowProto":
    constexpr static SpS OxidTankBtmProto =
      SpS
      (
        TT::UnDef(),
        false,                      // Facing Down
        OxidTankLowProto.GetLow()[0],
        MinD,
        OxidTankBtmH,
        Propellants::LOxDens
      );
    static_assert(OxidTankLowProto.GetLow()[0] == OxidTankBtmProto.GetUp()[0]);

    //-----------------------------------------------------------------------//
    // FuelTank Protos:                                                      //
    //-----------------------------------------------------------------------//
    // The Cylindrical Enclosure of the OxidTankBtm and FuelTankTop SpherSegms.
    // Its UpperBase is the Base of the "OxidTankBtmProto". It does not belong
    // to Oxid or Fuel Tank:
    constexpr static TrC OxidFuelEnclProto =
      TrC
      (
        TT::UnDef(),
        OxidTankBtmEnclProto.GetLow()[0],
        MinD,
        OxidFuelEnclH,              // No propellant there
        Density(0.0)
      );

    // The Top SpherSegm, with the Base coinciding with the LowBase of the
    // "OxidFuelEnclProto":
    constexpr static SpS FuelTankTopProto =
      SpS
      (
        TT::UnDef(),
        true,                       // Facing Up!
        OxidFuelEnclProto.GetLow()[0],
        MinD,
        FuelTankTopH,
        Propellants::RG1Dens
      );
    static_assert((OxidTankBtmProto.GetLow()[0] - FuelTankTopProto.GetUp()[0])
                  .ApproxEquals(OxidFuelGap));

    // The Mid Cylinder:
    constexpr static TrC FuelTankMidProto =
      TrC
      (
        TT::UnDef(),
        OxidFuelEnclProto.GetLow()[0],
        MinD,
        FuelTankMidH,
        Propellants::RG1Dens
      );

    // The Bottom SpherSegm, with the Base coinciding with the LowBase of the
    // "FuelTankMidProto":
    constexpr static SpS FuelTankBtmProto =
      SpS
      (
        TT::UnDef(),
        false,                      // Facing Down
        FuelTankMidProto.GetLow()[0],
        MinD,
        FuelTankBtmH,
        Propellants::RG1Dens
      );

    //-----------------------------------------------------------------------//
    // H2O2 Tank Proto:                                                      //
    //-----------------------------------------------------------------------//
    // The LowCylinder part extending beyond the FuelTankMid. It serves as the
    // Enclosure for FuelTankBtm, H2O2Tank, LiqN2Tank and the Engine:
    //
    constexpr static TrC TailEnclProto =
      TrC
      (
        TT::UnDef(),
        FuelTankMidProto.GetLow()[0],
        MinD,
        LowCylH - (OxidTankBtmEnclProto.GetLow()[0] -
                   FuelTankMidProto    .GetLow()[0]),
        Density(0.0)                // No propellant there
      );

    // The Top Half-Torical Segment:
    constexpr static Tor  H2O2TankTopProto =
      Tor
      (
        TT::UnDef(),
        true,                       // Facing Up
        FuelTankBtmProto.GetLow()[0] - H2O2Gap - H2O2TankMinoR,
        2.0 * H2O2TankMinoR,
        MinD,
        Propellants::H2O2Dens
      );

    // The DoubleCylinder Mid:
    constexpr static DCyl H2O2TankMidProto =
      DCyl
      (
        TT::UnDef(),
        H2O2TankTopProto.GetLow()[0],
        MinD,
        MinD - 2.0 * H2O2TankMinoR,
        H2O2TankCylH,
        Propellants::H2O2Dens
      );

    // The Bottom Half-Torical Segment:
    constexpr static Tor  H2O2TankBtmProto =
      Tor
      (
        TT::UnDef(),
        false,                      // Facing Down
        H2O2TankMidProto.GetLow()[0],
        2.0 * H2O2TankMinoR,
        MinD,
        Propellants::H2O2Dens
      );

    //-----------------------------------------------------------------------//
    // N2 Tank Proto:                                                        //
    //-----------------------------------------------------------------------//
    // Just 2 Semi-Toric Segments:
    //
    constexpr static Tor LiqN2TankTopProto =
      Tor
      (
        TT::UnDef(),
        true,                       // Facing Up
        H2O2TankBtmProto.GetLow()[0] - H2O2Gap - LiqN2TankMinoR,
        2.0 * LiqN2TankMinoR,
        MinD,
        Propellants::LN2Dens
      );

    constexpr static Tor LiqN2TankBtmProto =
      Tor
      (
        TT::UnDef(),
        false,                      // Facing Down
        LiqN2TankTopProto.GetLow()[0],
        2.0 * LiqN2TankMinoR,
        MinD,
        Propellants::LN2Dens
      );

    //-----------------------------------------------------------------------//
    // Scale Factor to determine the individual masses:                      //
    //-----------------------------------------------------------------------//
    // XXX: Once again, InterStageGrid is not modeled, its mass is spread over
    // the masses of other "ME"s:
    //
    constexpr static double ScaleFactor =
      ME::GetMassScale
      (
        // Components for which we provide the TotalMass (below):
        { &DeflectorConeProto,    &EquipBayProto,     &OxidTankTopEnclProto,
          &OxidTankTopProto,      &OxidTankUpProto,   &OxidTankLowProto,
          &OxidTankBtmProto,
          &OxidTankBtmEnclProto,  &OxidFuelEnclProto,
          &FuelTankTopProto,      &FuelTankMidProto,  &FuelTankBtmProto,
          &TailEnclProto,
          &H2O2TankTopProto,      &H2O2TankMidProto,  &H2O2TankBtmProto,
          &LiqN2TankTopProto,     &LiqN2TankBtmProto
        },
        // The TotalMass of the sbove MEs (unlike "Stage3", it does NOT include
        // Gases):
        EmptyMass - EngMass
      );

  public:
    //=======================================================================//
    // "MechElement"s with Real Masses:                                      //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Top of Stage2:                                                        //
    //-----------------------------------------------------------------------//
    constexpr static TrC DeflectorCone    =
      ME::ProRateMass(DeflectorConeProto,   ScaleFactor);

    constexpr static TrC EquipBay         =
      ME::ProRateMass(EquipBayProto,        ScaleFactor);

    //-----------------------------------------------------------------------//
    // OxidTank-Related :                                                    //
    //-----------------------------------------------------------------------//
    constexpr static TrC OxidTankTopEncl  =
      ME::ProRateMass(OxidTankTopEnclProto, ScaleFactor);

    constexpr static SpS OxidTankTop      =
      ME::ProRateMass(OxidTankTopProto,     ScaleFactor);

    constexpr static TrC OxidTankUp       =
      ME::ProRateMass(OxidTankUpProto,      ScaleFactor);

    constexpr static TrC OxidTankLow      =
      ME::ProRateMass(OxidTankLowProto,     ScaleFactor);

    constexpr static SpS OxidTankBtm      =
      ME::ProRateMass(OxidTankBtmProto,     ScaleFactor);

    constexpr static TrC OxidTankBtmEncl  =
      ME::ProRateMass(OxidTankBtmEnclProto, ScaleFactor);

    //-----------------------------------------------------------------------//
    // FuelTank-Related (incl the OxidFuelEncl):                             //
    //-----------------------------------------------------------------------//
    constexpr static TrC OxidFuelEncl     =
      ME::ProRateMass(OxidFuelEnclProto,    ScaleFactor);

    constexpr static SpS FuelTankTop      =
      ME::ProRateMass(FuelTankTopProto,     ScaleFactor);

    constexpr static TrC FuelTankMid      =
      ME::ProRateMass(FuelTankMidProto,     ScaleFactor);

    constexpr static SpS FuelTankBtm      =
      ME::ProRateMass(FuelTankBtmProto,     ScaleFactor);

    //-----------------------------------------------------------------------//
    // Tail Section:                                                         //
    //-----------------------------------------------------------------------//
    constexpr static TrC  TailEncl        =
      ME::ProRateMass(TailEnclProto,       ScaleFactor);

    constexpr static Tor  H2O2TankTop     =
      ME::ProRateMass(H2O2TankTopProto,    ScaleFactor);

    constexpr static DCyl H2O2TankMid     =
      ME::ProRateMass(H2O2TankMidProto,    ScaleFactor);

    constexpr static Tor  H2O2TankBtm     =
      ME::ProRateMass(H2O2TankBtmProto,    ScaleFactor);

    constexpr static Tor  LiqN2TankTop    =
      ME::ProRateMass(LiqN2TankTopProto,   ScaleFactor);

    constexpr static Tor  LiqN2TankBtm    =
      ME::ProRateMass(LiqN2TankBtmProto,   ScaleFactor);

    // The following invariant must hold:
    static_assert
      ((LiqN2TankBtm.GetLow()[0] - EngineH).ApproxEquals(NozzlesLowX));

    // The lowest point of the TailEncl (still not of the whole Stage2, because
    // Engine Nozzles extend beyond the TailEncl):
    constexpr static Len TailEnclLowX     = TailEncl.GetLow()[0];
    static_assert(NozzlesLowX < TailEnclLowX);

    // RD-108A Engine, modeled as a PointMass.   XXX: We assume that its CoM is
    // located at 1/3 height from the Engine Top, but this is obviously GROSSLY
    // IMPRECISE:
    constexpr static Len EngineCoMX       =
      LiqN2TankBtm.GetLow()[0] - EngineH / 3.0;

    constexpr static PM  Engine =
      PM
      (
        TT::UnDef(),
        EngineCoMX,
        0.0_m,
        0.0_m,
        EngMass
      );

    //-----------------------------------------------------------------------//
    // Empty Stage2 as a whole:                                              //
    //-----------------------------------------------------------------------//
    // Unlike Stage3, we do NOT include the Pressurisation Gas, in this case
    // N2, into the "EmptyME":
    //
    constexpr static ME EmptyME      =
      DeflectorCone   + EquipBay     + OxidTankTopEncl +
      OxidTankTop     + OxidTankUp   + OxidTankLow     + OxidTankBtm +
      OxidTankBtmEncl + OxidFuelEncl +
      FuelTankTop     + FuelTankMid  + FuelTankBtm     +
      TailEncl        +
      H2O2TankTop     + H2O2TankMid  + H2O2TankBtm     +
      LiqN2TankTop    + LiqN2TankBtm +
      Engine;

    // The over-all EmptyMass check:
    static_assert(EmptyME.GetMass().ApproxEquals(EmptyMass));

  private:
    //=======================================================================//
    // Propellant Volumes and Mass Capacities:                               //
    //=======================================================================//
    // Propallant Mass Capacities (MC) of Oxid, Fuel, H2O2 and N2 Tank Sections
    // and their "Unions".
    // Oxid:
    constexpr static Mass OxidTankTopMC      = OxidTankTop.GetPropMassCap();
    constexpr static Mass OxidTankUpMC       = OxidTankUp .GetPropMassCap();
    constexpr static Mass OxidTankLowMC      = OxidTankLow.GetPropMassCap();
    constexpr static Mass OxidTankBtmMC      = OxidTankBtm.GetPropMassCap();
    // Unions:
    constexpr static Mass OxidTankBtmLowMC   = OxidTankBtmMC    + OxidTankLowMC;
    constexpr static Mass OxidTankBtmLowUpMC = OxidTankBtmLowMC + OxidTankUpMC;

    // Fuel:
    constexpr static Mass FuelTankTopMC      = FuelTankTop.GetPropMassCap();
    constexpr static Mass FuelTankMidMC      = FuelTankMid.GetPropMassCap();
    constexpr static Mass FuelTankBtmMC      = FuelTankBtm.GetPropMassCap();
    // Union:
    constexpr static Mass FuelTankBtmMidMC   = FuelTankBtmMC    + FuelTankMidMC;

    // H2O2:
    constexpr static Mass H2O2TankTopMC      = H2O2TankTop.GetPropMassCap();
    constexpr static Mass H2O2TankMidMC      = H2O2TankMid.GetPropMassCap();
    constexpr static Mass H2O2TankBtmMC      = H2O2TankBtm.GetPropMassCap();
    // Union:
    constexpr static Mass H2O2TankBtmMidMC   = H2O2TankBtmMC    + H2O2TankMidMC;

    // LiqN2:
    constexpr static Mass LiqN2TankTopMC     = LiqN2TankTop.GetPropMassCap();
    constexpr static Mass LiqN2TankBtmMC     = LiqN2TankBtm.GetPropMassCap();

  public:
    // Maximum Theoretical Capacities of the resp Tanks:
    constexpr static Mass OxidTankMC  = OxidTankBtmLowUpMC + OxidTankTopMC;
    constexpr static Mass FuelTankMC  = FuelTankBtmMidMC   + FuelTankTopMC;
    constexpr static Mass H2O2TankMC  = H2O2TankBtmMidMC   + H2O2TankTopMC;
    constexpr static Mass LiqN2TankMC = LiqN2TankTopMC     + LiqN2TankBtmMC;

    // Propellant Load Ratios (ActualMass / TheorMassCapacity):
    constexpr static double OxidLoadRatio  = double(OxidMass  / OxidTankMC);
    static_assert(OxidLoadRatio  < 1.0);
    constexpr static double FuelLoadRatio  = double(FuelMass  / FuelTankMC);
    static_assert(FuelLoadRatio  < 1.0);
    constexpr static double H2O2LoadRatio  = double(H2O2Mass  / H2O2TankMC);
    static_assert(H2O2LoadRatio  < 1.0);
    constexpr static double LiqN2LoadRatio = double(LiqN2Mass / LiqN2TankMC);
    static_assert(LiqN2LoadRatio < 1.0);

  private:
    //-----------------------------------------------------------------------//
    // "ME" objs for Max Theoretical Propellant Loads in Tank Sections:      //
    //-----------------------------------------------------------------------//
    // (For optimisation of "GetDynParams"):
    // Oxid:
    constexpr static ME OxidTopME         = OxidTankTop.GetPropBulkME();
    constexpr static ME OxidUpME          = OxidTankUp .GetPropBulkME();
    constexpr static ME OxidLowME         = OxidTankLow.GetPropBulkME();
    constexpr static ME OxidBtmME         = OxidTankBtm.GetPropBulkME();
    // Unions:
    constexpr static ME OxidBtmLowME      = OxidBtmME      + OxidLowME;
    constexpr static ME OxidBtmLowUpME    = OxidBtmLowME   + OxidUpME;
    // The whole maximum Oxid bulk:
    constexpr static ME OxidME            = OxidBtmLowUpME + OxidTopME;

    //Fuel:
    constexpr static ME FuelTopME         = FuelTankTop.GetPropBulkME();
    constexpr static ME FuelMidME         = FuelTankMid.GetPropBulkME();
    constexpr static ME FuelBtmME         = FuelTankBtm.GetPropBulkME();
    // Union:
    constexpr static ME FuelBtmMidME      = FuelBtmME      + FuelMidME;
    // The whole maximum Fuel bulk:
    constexpr static ME FuelME            = FuelBtmMidME   + FuelTopME;

    // For H2O2 and LiqN2, "Top" MEs are not required:
    // H2O2:
    constexpr static ME H2O2TopME         = H2O2TankTop.GetPropBulkME();
    constexpr static ME H2O2MidME         = H2O2TankMid.GetPropBulkME();
    constexpr static ME H2O2BtmME         = H2O2TankBtm.GetPropBulkME();
    // Union:
    constexpr static ME H2O2BtmMidME      = H2O2BtmME      + H2O2MidME;
    // The whole maximum H2O2 bulk:
    constexpr static ME H2O2ME            = H2O2BtmMidME   + H2O2TopME;

    // LiqN2:
    constexpr static ME LiqN2TopME        = LiqN2TankTop.GetPropBulkME();
    constexpr static ME LiqN2BtmME        = LiqN2TankBtm.GetPropBulkME();
    // The whole maximum H2O2 bulk:
    constexpr static ME LiqN2ME           = LiqN2TopME     + LiqN2BtmME;

  public:
    //=======================================================================//
    // Thrust Vector Control:                                                //
    //=======================================================================//
    // Achieved via 4 Gimbaled Vernier Chambers. Each Vernier Chamber is located
    // in the XY or XZ plane and gimbaled around the corresp axis:
    //   0: @ -Y (towards Block B, Angle > 0 ==> towards -Z),
    //   1: @ -Z (towards Block V, Angle > 0 ==> towards +Y),
    //   2: @ +Y (towards Block G, Angle > 0 ==> towards +Z),
    //   3: @ +Z (towards Block D, Angle > 0 ==> towards -Y),
    // and gimbaling results in "positive" thrust projections in the opposite
    // directions (+Z, -Y, -Z, +Y, resp).
    // Each Vernier can be deflected within the GimbalAmpl in the Tangential
    // Plane. Similar to Stage3, we assume that Gimbaling (Deflection) Angle
    // is positive when the corresp Vernier is deflected  Counter-Clock-Wise
    // in the YZ plane, viwed in the right-handed XYZ co-ord system, ie from
    // the positive end of the OX axis. And again, similar to Stage3, opposite
    // Verniers would normally be moved symmetrically in one direction, so that
    //    VernDeflections[0] = - VernDeflections[2],
    //    VernDeflections[1] = - VernDeflections[3],
    // but this is not strictly enforced:
    //
    constexpr static Angle_deg VernGimbalAmpl = Angle_deg(45.0);
    using            VernDeflections          = std::array<Angle_deg, 4>;

    //=======================================================================//
    // Dynamic Params as a function of Flight Time:                          //
    //=======================================================================//
    // NB: This method is NOT "constexpr": it is intended to be called at Run-
    // Time (eg multiple times during Trajectory Integration):
    //
    static StageDynParams<LVSC::Soyuz21b>
    GetDynParams
    (
      FT                     a_ft,
      Pressure               a_p,      // Curr Atmospheric Pressure
      VernDeflections const& a_vern_defls
    );
  };
}
// End namespace SpaceBallistics
