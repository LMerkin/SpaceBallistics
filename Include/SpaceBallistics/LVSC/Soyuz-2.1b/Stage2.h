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
#include <cassert>

namespace SpaceBallistics
{
  namespace SC  = Soyuz21b_Consts;
  namespace SH  = Soyuz21b_Head;

  // All "MechElements" are instantiated with "LVSC::Soyuz21b":
  using     ME   = MechElement   <LVSC::Soyuz21b>;
  using     TrC  = TrCone        <LVSC::Soyuz21b>;
  using     SpS  = SpherSegm     <LVSC::Soyuz21b>;
  using     Tor  = ToricSegm     <LVSC::Soyuz21b>;
  using     DCyl = DoubleCylinder<LVSC::Soyuz21b>;

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
    constexpr static Len    TopX          = SC::X0 - SC::Stage3Len;

    // The IntegrStageGrid:
    constexpr static Len    GridH         = 0.975_m;

    // DeflectorCone Height and Base Diameter. The Height is < GridH, so the
    // DeflectorCone is within the Grid:
    constexpr static Len    DeflConeH     = GridH / 2.0;
    constexpr static Len    TopD          = 2.700_m;  // Slightly >than Stage3D

    // Equipment Bay: Cylindrical or expanding TrCone? Information differs in
    // different drawings. Assume Cylindrical, of the same "TopD" diameter:
    constexpr static Len    EquipBayH     = 0.595_m;

    //-----------------------------------------------------------------------//
    // OxidTank:                                                             //
    //-----------------------------------------------------------------------//
    // The over-all Upper (Expanding) TrCone: Contains the Top SpherSegm of the
    // OxidTank, and the Upper TrCone of the OxidTank itself. The diameter inc-
    // reases from "TopD" to "MaxD":
    constexpr static Len    UpperTrCH     = 5.840_m;  // Or 5.835 m?
    constexpr static Len    MaxD          = 2.950_m;  // Max for entire Stage2

    // The over-all Middle (Narrowing-Down) TrCone: Contains the Low TrCone of
    // the OxidTank. The diameter decreses from MaxD to MinD:
    constexpr static Len    MidTrCH       = 4.55_m;
    constexpr static Len    MinD          = 2.05_m;

    // Top OxidTank SpherSegm. Its top-most point (the pole) is at the level of
    // the notional bottom base of the EquipBay:
    constexpr static Len    OxidTankTopH  = 0.8_m;    // Approx
    constexpr static Len    OxidTankTopD  =
      TopD + double(OxidTankTopH / UpperTrCH) * (MaxD - TopD);

    // Height of the Upper TrCone part of the OxidTank: extends to the end of
    // the over-all  Upper TrCone:
    constexpr static Len    OxidTankUpH   = UpperTrCH - OxidTankTopH;

    // Height of the Lower TrCone part of the OxidTank (extends ALMOST to the
    // end of the Middle TrCone, so its diameter is slightly > MinD):
    constexpr static Len    OxidTankBtmH  = 0.5_m;    // Approx
    constexpr static Len    OxidTankBtmH2 = OxidTankBtmH / 2.0;
    constexpr static Len    OxidTankLowH  = MidTrCH   - OxidTankBtmH2;
    constexpr static Len    OxidTankLowD  =
      MinD + (MaxD - MinD) * double(OxidTankBtmH2 / MidTrCH);

    // Bottom SpherSegm of OxidTank: D=OxidTankLowD, H=OxidTankBtmH

    //-----------------------------------------------------------------------//
    // FuelTank:                                                             //
    //-----------------------------------------------------------------------//
    // Located entirely in the Cylindrical Low section of Stage2. The Diameter
    // of this section is "MinD", and over-all length is:
    constexpr static Len    LowCylH       = 15.178_m;

    // There is a gap between the OxidTankBtm and FuelTankTop. Its exact size
    // is not known, so assume something similar to Stage3:
    constexpr static Len    OxidFuelGap   = 0.15_m;

    // Top SpherSegm:
    constexpr static Len    FuelTankTopH  = OxidTankBtmH;   // XXX: Approx...

    // Bottom SpherSegm:
    constexpr static Len    FuelTankBtmH  = FuelTankTopH;

    // TrCone Enclosure of the top-half of "OxidTankBtm" (NB: the lower-half is
    // enclosed by the "OxidFuelEnc", see below):
    constexpr static Len    OxidBtmEnclH  = OxidTankBtmH2;

    // Length of the main cylindrical section of the FuelTank: Known only
    // approximately. We assume (see above) that half of the OxidTankBtm secti-
    // on is in the "MidTrCone", and the "FuelTankBtm" SpherSegm is beyond the
    // 10.378_m quoted below:
    constexpr static Len    OxidFuelEnclH =
      OxidTankBtmH2 + OxidFuelGap + FuelTankTopH;

    constexpr static Len    FuelTankMidH  = 10.378_m - OxidFuelEnclH;

    //-----------------------------------------------------------------------//
    // H2O2 and N2 Tanks:                                                    //
    //-----------------------------------------------------------------------//
    // Top and Bottom parts are Semi-ToricSegms, the Mid part is a Double
    // Cylinder. All sizes are approximate:
    constexpr static Len    H2O2TankOutR  = MinD / 2.0;
    constexpr static Len    H2O2TankInR   = 0.63_m;
    constexpr static Len    H2O2TankMinoR = (H2O2TankOutR - H2O2TankInR) / 2.0;
    constexpr static Len    H2O2TankCylH  = 1.2_m;

    // N2 Tank is a full Torus:
    constexpr static Len    N2TankMinoR   = 0.2_m;

    // We assume that the two Gaps (FuelTank-H2O2Tank and H2O2Tank -- N2Tank)
    // are the same, and derive them from the over-all compartment length. In
    // that part of the Tail compartment, the following are located:
    //   FuelTankBtm, H2O2Gap, H2O2Tank, H2O2Gap, N2Tank:
    // Make this gap pretty small, otherwise the Engine may not fit:
    constexpr static Len    H2O2Gap       = 0.05_m;

    //-----------------------------------------------------------------------//
    // Over-all length of Stage2:                                            //
    //-----------------------------------------------------------------------//
    // (From the top of the InterStageGrid to the end of the LowCyl, but w/o
    // the extending Engine Nozzles):
    constexpr static Len    H =
      GridH + EquipBayH + UpperTrCH + MidTrCH + LowCylH;
    static_assert(H.ApproxEquals(27.138_m));

    // Nozzles extension beyond "H":
    constexpr static Len   EngineNozzlesExtH  = 0.635_m;

    // RD-108A Height:
    constexpr static Len   EngineH            = 2.865_m;

    //=======================================================================//
    // Masses:                                                               //
    //=======================================================================//
    // EmptyMass: XXX: StarSem says 6545 kg:
    constexpr static Mass   EmptyMass     = 6450.0_kg;

    // Mass of the RD-108A engine (part of EmptyMass):
    constexpr static Mass   EngMass       = 1075.0_kg;

    // Masses of Fuel (Naftil) and Oxidiser. As for Stage3, FuelMass includes
    // extra 0.2% for the antifreeze (2-EtoxyEthanol):
    constexpr static Mass   FuelMass      = 26794.0_kg * 1.002;
                                                        // StarSem: 26300 (T1)
    constexpr static Mass   OxidMass      = 63709.0_kg; // StarSem: 63800
    constexpr static Mass   H2O2Mass      = 2636.0_kg;
    constexpr static Mass   N2Mass        = 513.0_kg;
    constexpr static Mass   FullMass      = EmptyMass + FuelMass + OxidMass +
                                            H2O2Mass  + N2Mass;

    // UnSpendable Remnants of the Fuel and Oxidiser in Stage2 at the engine
    // cut-off time.   NB: The Remnants are about 1% of the corresp initial
    // masses. NB: They are Technically UnSpendable, so they do NOT include the
    // "guarantee margins":
    constexpr static Mass   FuelRem       = 272.0_kg;
    constexpr static Mass   OxidRem       = 678.0_kg;
    constexpr static Mass   H2O2Rem       = 263.0_kg;
    constexpr static Mass   N2Rem         = 88.0_kg;

    //-----------------------------------------------------------------------//
    // RD-108A (14D21) Engine Performance:                                   //
    //-----------------------------------------------------------------------//
    // Isp (SL/Vac, sec): EnergoMash says 257.7/320.6,  StarSem: 255.0/319.0;
    // but it is unclear whether these vals apply to the Main Engine    (w/o 
    // Vernier Chambers) or to the Whole Engine (incl Vernier Chambers, which
    // may have a slightly lower Isp). So assume higher vals for the Main Eng-
    // ine:
    constexpr static Time   IspMainSL     = 262.9_sec;
    constexpr static Time   IspMainVac    = 327.1_sec;

    // XXX: For Thrust (SL/Vac, tf), many different vals exist:
    // 70/87, 79.2/92.1, 80.8/94 (EnergoMash), 80.81/100.97;
    // but we must have IspVac/IspSL = ThrustVac/ThrustSL = 1.244  for the Main
    // Engine (see above);   then the 2nd and 3rd vals are clearly incorrect
    // (the ratio for them is 1.163,  so they are probably made of a mix of the
    // Main Engine and Vernier Chamber vals); the 1st and 4th vals yield reas-
    // onable ratios;  then the 1st one is apparently for Main Engine only,
    // and the last one is for the engine over-all;   use the former here:
    //
    constexpr static Force  ThrustMainSL  = 70000.0_kg * g0;
    constexpr static Force  ThrustMainVac =
      ThrustMainSL * double(IspMainVac / IspMainSL);  // 87.1 tf

    // For each Vernier Chamber, the data are taken from the old RD-107, so
    // they might not be very accurate. For example, the RD-107 Vernier Chamber
    // Fuel/Oxid Rates of 4.15/8.55 kg/sec would then lead to thrust of 3.2/4.0
    // tf, and the total Engine thrust (Main + 4 Vernier Chambers) of  82.8/103
    // tf, which is higher by 2 tf than the above figure of 80.81/100.97; so we
    // keep the Isp vals (their proportion is almost the same as for the  Main
    // Engine) and adjust the MassRate:
    // IMPORTANT:
    // Vernier Chambers of Stage2 are installed in positions at +-Y, +-Z  axes,
    // ie towards the corresp Stage1 blocks. They are rotatable to +-45 degress
    // in the planes parallel to the X axis, ie, the Vernier Chambers installed
    // at +-Y are rotatable in the XZ plane, and those installed at +-Z are ro-
    // tatable in the XY plane:
    //
    constexpr static Time   IspVernSL1     = 251.9_sec;
    constexpr static Time   IspVernVac1    = 313.1_sec;
    constexpr static Force  ThrustVernSL1  = 2700.0_kg * g0;
    constexpr static Force  ThrustVernVac1 =
      ThrustVernSL1 * double(IspVernVac1 / IspVernSL1);
      // 3.356 tf => Total=80.8/100.5 tf, very close to 80.81/100.97 tf above!

    // We can then calculate the MassRates for the Main Engine and for each
    // Vernier Chamber:
    constexpr static MassRate MassRateMain = ThrustMainSL / (IspMainSL  * g0);
    static_assert(MassRateMain.ApproxEquals(ThrustMainVac / (IspMainVac * g0)));
    // ~277.89 kg/sec

    constexpr static MassRate MassRateVern1 =
      ThrustVernSL1 / (IspVernSL1 * g0);
    static_assert
      (MassRateVern1.ApproxEquals(ThrustVernVac1 / (IspVernVac1 * g0)));
    //
    // ~10.72  kg/sec, NOT 4.15+8.55=12.70 kg/sec as it would be with the orig-
    // inal RD-107 Vernier Chamber data...

    constexpr static MassRate MassRateVern4 = 4.0 * MassRateVern1;

    // Total MassRate for the whole Engine:
    constexpr static MassRate MassRateEng  = MassRateMain + MassRateVern4;
    // ~320.76 kg/sec

    // Separate Fuel and Oxid Rates are obtained using the Oxid/Fuel Ratio which
    // we derive from the over-all Fuel and Oxid spendable masses:
    constexpr static double OxidFuelRatio =
      double((OxidMass - OxidRem) / (FuelMass - FuelRem));
    // ~2.38, very close to 2.39 often quoted for RD-108A

    constexpr static MassRate MassRateFuel =
      MassRateEng                  / (1.0  + OxidFuelRatio);
    constexpr static MassRate MassRateOxid =
      MassRateEng * (OxidFuelRatio / (1.0  + OxidFuelRatio));

    // IMPORTANT: For RD-108A, we assume that FullThrust instant is the same as
    // LiftOff (Contact Separation), ie t=0.
    // Preliminary thrust level (assumed to be a 25%) occurs notionally at -15
    // sec (in reality, the ignition sequence is more complex):
    // So the Propellant Mass spent at Ignition is:
    // XXX: It is OK to use hard-wired consts here, they are not used anywhere
    // else (unlike the ShutDown process, see below):
    //
    constexpr static Mass     IgnPropMass = MassRateEng * 15.0_sec * 0.25;

    // RD-108A Shut-Down Sequence:
    // If "tF" is the Full ShutDown Time,  then the Main Chambers run for 1 sec
    // (from tF-6 to tF-5) sec at the 25% (???) thrust level and then shut down
    // completely, and the Vernier Chambers run (at which thrust level? 25% as
    // well?) until "tF".
    // Thus, the Propellany Mass spent at ShutDown is:
    //
    constexpr static Time     VernThrottlAdvance = 6.0_sec;
    constexpr static Time     MainThrottlAdvance = 1.0_sec;
    constexpr static double   ThrottlLevel       = 0.25;
    constexpr static Mass     ShutDownPropMass   =
      (MassRateMain * MainThrottlAdvance + MassRateVern4 * VernThrottlAdvance) *
      ThrottlLevel;

    // Thus, the Propellant Mass left for the FullThrust Mode:
    constexpr static Mass     FullThrustPropMass =
      FuelMass + OxidMass - IgnPropMass - ShutDownPropMass - FuelRem - OxidRem;
    static_assert(IsPos(FullThrustPropMass));

    // Then the Max Time at FullThrust is:
    constexpr static Time     MaxFullThrustTime  =
      FullThrustPropMass / MassRateEng;

    // Then the Max Time of RD-108A operation from LiftOff=FullThrust (NOT from
    // Ignition which occurs before LiftOff) to Full ShutDown  is.
    // XXX: For Stage3, a similar param is called "MaxBurnDur", but in this case
    // the ignition occurs BEFORE the lift-off, and we are actually interested
    // in the "MaxFlightTime", not "MaxBurnDur":
    constexpr static Time     MaxFlightTime      =
      MaxFullThrustTime + std::max(VernThrottlAdvance, MainThrottlAdvance);

    // MaxFlightTime appears to be ~291 sec, whereas the actual Stage2CutOffTime
    // is ~286 sec, which is a reasonably good match:
    static_assert(SC::Stage2CutOffTime < MaxFlightTime);

    // NB: In addition, there is a MassRate due to H2O2 burning   to drive the
    // TurboPumps; however, this MassRate does not formally participate in the
    // Thrust. We assume that H2O2 flow begins at FullThrust / LiftOff time
    // and lasts for the whole "MaxFlightTime":
    //
    constexpr static MassRate MassRateH2O2  =
      (H2O2Mass - H2O2Rem) / MaxFlightTime;

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
    // Enclosure for FuelTankBtm, H2O2Tank, N2Tank and the Engine:
    //
    constexpr static TrC TailEnclProto =
      TrC
      (
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
    constexpr static Tor N2TankTopProto =
      Tor
      (
        true,                       // Facing Up
        H2O2TankBtmProto.GetLow()[0] - H2O2Gap - N2TankMinoR,
        2.0 * N2TankMinoR,
        MinD,
        Propellants::LN2Dens
      );
    constexpr static Tor N2TankBtmProto =
      Tor
      (
        false,                      // Facing Down
        N2TankTopProto.GetLow()[0],
        2.0 * N2TankMinoR,
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
          &N2TankTopProto,        &N2TankBtmProto
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

    constexpr static Tor  N2TankTop       =
      ME::ProRateMass(N2TankTopProto,      ScaleFactor);

    constexpr static Tor  N2TankBtm       =
      ME::ProRateMass(N2TankBtmProto,      ScaleFactor);

    // Two ways of computing the Lowest Co-Ord (the Nozzles End), should be
    // equal up to 3 cm:
    constexpr static Len  EngineNozzlesLow1 =
      N2TankBtm.GetLow()[0] - EngineH;
    constexpr static Len  EngineNozzlesLow2 =
      TopX - H - EngineNozzlesExtH;
    static_assert(Abs(EngineNozzlesLow1 - EngineNozzlesLow2) < 0.03_m);

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
    constexpr static Mass H2O2TankTopMC     = H2O2TankTop.GetPropMassCap();
    constexpr static Mass H2O2TankMidMC     = H2O2TankMid.GetPropMassCap();
    constexpr static Mass H2O2TankBtmMC     = H2O2TankBtm.GetPropMassCap();
    // Union:
    constexpr static Mass H2O2TankBtmMidMC  = H2O2TankBtmMC    + H2O2TankMidMC;

    // N2:
    constexpr static Mass N2TankTopMC       = N2TankTop.GetPropMassCap();
    constexpr static Mass N2TankBtmMC       = N2TankBtm.GetPropMassCap();

  public:
    // Volumes of the Oxid, Fuel, H2O2 and N2 Tanks:
    constexpr static Vol  OxidTankVol  =
      OxidTankTop.GetEnclVol() + OxidTankUp .GetEnclVol() +
      OxidTankLow.GetEnclVol() + OxidTankBtm.GetEnclVol();

    // Maximum Theoretical Capacities of the resp Tanks:
    constexpr static Mass OxidTankMC   = OxidTankBtmLowUpMC + OxidTankTopMC;
    constexpr static Mass FuelTankMC   = FuelTankBtmMidMC   + FuelTankTopMC;
    constexpr static Mass H2O2TankMC   = H2O2TankBtmMidMC   + H2O2TankTopMC;
    constexpr static Mass N2TankMC     = N2TankTopMC        + N2TankBtmMC;

    // Propellant Load Ratios (ActualMass / TheorMassCapacity):
    constexpr static double OxidLoadRatio = double(OxidMass / OxidTankMC);
    static_assert(OxidLoadRatio < 1.0);
    constexpr static double FuelLoadRatio = double(FuelMass / FuelTankMC);
    static_assert(FuelLoadRatio < 1.0);
    constexpr static double H2O2LoadRatio = double(H2O2Mass / H2O2TankMC);
    static_assert(H2O2LoadRatio < 1.0);
    constexpr static double N2LoadRatio   = double(N2Mass   / N2TankMC);
    static_assert(N2LoadRatio   < 1.0);

    //=======================================================================//
    // Thrust Vector Control:                                                //
    //=======================================================================//
    // Achieved via 4 Gimbaled Vernier Chambers. Each Vernier Chamber is located
    // in the XY or XZ plane:
    //   0: @ -Y
    //   1: @ -Z
    //   2: @ +Y
    //   3: @ +Z
    // Each chamber can be deflected within the GimbalAmpl in the Tangential
    // Plane. Similar to Stage3, we assume that Gimbaling (Deflection) Angle
    // is positive when the corresp chamber is deflected  Counter-Clock-Wise
    // in the YZ plane. And again, similar to Stage3, opposite Chambers would
    // normally be moved symmetrically in one direction, so that
    //    ChamberDeflections[0] = - ChamberDeflections[2],
    //    ChamberDeflections[1] = - ChamberDeflections[3],
    // but this is not strictly enforced:
    //
    constexpr static Angle_deg GimbalAmpl = Angle_deg(45.0);
    using            ChamberDeflections   = std::array<Angle_deg, 4>;

    //=======================================================================//
    // Dynamic Params as a function of Flight Time:                          //
    //=======================================================================//
    // NB: This method is NOT "constexpr": it is intended to be called at Run-
    // Time (eg multiple times during Trajectory Integration):
    //
    static StageDynParams<LVSC::Soyuz21b>
    GetDynParams(Time a_t, ChamberDeflections const& a_chamber_defls);
  };
}
// End namespace SpaceBallistics
