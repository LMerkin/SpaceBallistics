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
  using     ME  = MechElement<LVSC::Soyuz21b>;
  using     TC  = TrCone     <LVSC::Soyuz21b>;
  using     SpS = SpherSegm  <LVSC::Soyuz21b>;

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
    // Consts:                                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Geometry:                                                             //
    //-----------------------------------------------------------------------//
    // Stage2 Top (Base of Deflector Cone): Bottom of Stage3 minus the Height
    // of InterStage Grid:
    constexpr static Len    GridH         = 0.975_m;
    constexpr static Len    TopX          = SC::X0 - SC::Stage3Len - GridH;

    // DeflectorCone Height and Base Diameter. The Height is < GridH, so the
    // DeflectorCone is within the Grid:
    constexpr static Len    DeflConeH     = GridH / 2.0;
    constexpr static Len    TopD          = 2.700_m;  // Slightly >than Stage3D

    // Equipment Bay: Cylindrical or expanding TrCone? Information differs in
    // different drawings. Assume Cylindrical, of the same "TopD" diameter:
    constexpr static Len    EquipBayH     = 0.595_m;

    // The over-all Upper (Expanding) TrCone: Contains the Top SpherSegm of the
    // OxidTank, and the Upper TrCone of the OxidTank itself. The diameter inc-
    // reases from "TopD" to "MaxD":
    constexpr static Len    UpperTrCH     = 5.835_m;
    constexpr static Len    MaxD          = 2.950_m;  // Max for entire Stage2

    // The over-all Middle (Narrowing-Down) TrCone: Contains the Low TrCone and
    // Bottom SpherSegm of the OxidTank, as well as Top SpherSegm of the Fuel
    // Tank. The diameter decreses from MaxD to MinD:
    constexpr static Len    MidTrCH       = 4.55_m;
    constexpr static Len    MinD          = 2.05_m;

    // The Lower Cylinder: Contains the Mid Cylinder and the Bottom SpherSegm
    // of the FuelTank, the H2O2 Tank and the LiqN2 Tank (see below).

    // OxidTank:
    // Top SpherSegm. Its top-most point (the pole) is at the level of the not-
    // ional bottom base of the EqipBay:
    constexpr static Len    OxidTankTopH  = 0.900_m;  // Approx
    constexpr static Len    OxidTankTopD  =
      TopD + double(OxidTankTopH / UpperTrCH) * (MaxD - TopD);

    // Height of the Upper TrCone part of the OxidTank: Stretches to the end of
    // the over-all  Upper TrCone:
    constexpr static Len    OxidTankUpH   = UpperTrCH - OxidTankTopH;

    // Height of the Lower TrCone part of the OxidTank (located within the over-
    // all Middle TrCone):
    constexpr static Len    OxidTankLowH  = 3.20_m;

    // Height and Base Diameter of the Bottom SpherSegm part of the OxidTank:
    constexpr static Len    OxidTankBtmH  = 0.635_m;
    constexpr static Len    OxidTankBtmD  =
      MaxD + double(OxidTankLowH / MidTrCH) * (MinD - MaxD);

    // FuelTank:
    // Top SpherSegm (Base Diameter = MinD):
    constexpr static Len    FuelTankTopH  = 0.5_m;

    // There is a gap between the OxidTankBtm and FuelTankTop:
    constexpr static Len    InterTankGap  =
      MidTrCH - OxidTankLowH - OxidTankBtmH - FuelTankTopH;
    static_assert(IsPos(InterTankGap));

    // Length of the main cylindrical section of the FuelTank:
//  constexpr static Len    FuelTankMidH  =

    //-----------------------------------------------------------------------//
    // Masses:                                                               //
    //-----------------------------------------------------------------------//
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
    // cut-off time.   NB: The Remnanats are about 1% of the corresp initial
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
    constexpr static Time   IspVernSL      = 251.9_sec;
    constexpr static Time   IspVernVac     = 313.1_sec;
    constexpr static Force  ThrustVernSL   = 2700.0_kg * g0;
    constexpr static Force  ThrustVernVac  =
      ThrustVernSL * double(IspVernVac / IspVernSL);
      // 3.356 tf => Total=80.8/100.5 tf, very close to 80.81/100.97 tf above!

    // We can then calculate the MassRates for the Main Engine and for each
    // Vernier Chamber:
    constexpr static MassRate MassRateMain = ThrustMainSL / (IspMainSL  * g0);
    static_assert(MassRateMain.ApproxEquals(ThrustMainVac / (IspMainVac * g0)));
    // ~277.89 kg/sec

    constexpr static MassRate MassRateVern = ThrustVernSL / (IspVernSL  * g0);
    static_assert(MassRateVern.ApproxEquals(ThrustVernVac / (IspVernVac * g0)));
    // ~10.72  kg/sec, NOT 4.15+8.55=12.70 kg/sec as it would be with the orig-
    // inal RD-107 Vernier Chamber data...

    // Total MassRate for the whole Engine:
    constexpr static MassRate MassRateEng  = MassRateMain + 4.0 * MassRateVern;
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
    //
    constexpr static Mass     IgnPropMass = MassRateEng * 15.0_sec * 0.25;

    // RD-108A Shut-Down Sequence:
    // If "tF" is the Full ShutDown Time,  then the Main Chambers run for 1 sec
    // (from tF-6 to tF-5) sec at the 25% (???) thrust level and then shut down
    // completely, and the Vernier Chambers run (at which thrust level? 25% as
    // well?) until "tF".
    // Thus, the Propellany Mass spent at ShutDown is:
    //
    constexpr static Mass     ShutDownPropMass   =
      (MassRateMain * 1.0_sec + 4.0 * MassRateVern * 6.0_sec) * 0.25;

    // Thus, the Propellant Mass left for the FullThrust Mode:
    constexpr static Mass     FullThrustPropMass =
      FuelMass + OxidMass - IgnPropMass - ShutDownPropMass - FuelRem - OxidRem;
    static_assert(IsPos(FullThrustPropMass));

    // Then the Max Time at FullThrust is:
    constexpr static Time     MaxFullThrustTime  =
      FullThrustPropMass / MassRateEng;

    // Then the Max Time of RD-108A operation from LiftOff=FullThrust (NOT from
    // Ignition which occurs before LiftOff) to Full ShutDown  is:
    constexpr static Time     MaxFlightTime      = MaxFullThrustTime + 6.0_sec;

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
    //-----------------------------------------------------------------------//
    // "MechElement"s: "Proto"s with Yet-UnKnown Masses:                     //
    //-----------------------------------------------------------------------//
    //-----------------------------------------------------------------------//
    // Top-Level ME Protos:                                                  //
    //-----------------------------------------------------------------------//
    // XXX: Inter-Stage Grid (between Stages 2 and 3)  is  not modeled at all.
    // Its mass is included into the EmptyMass of Stage2, and will be "spread"
    // over the masses of all MEs defined below.  The reason is that the Grid
    // cannot (as yet) be represented as any supported "ME".  This will have a
    // minor effect on MoIs accuracy...
    //
    // "DeflectorCone" at the top of Stage2. Its base is @ "TopX" co-ord, but
    // the Ctor requires the vertex point in this case:
    //
    constexpr static TC DeflectorConeProto =
      TC
      (
        TopX - DeflConeH,               // Cone vertex point
        0.0_m,                          //  ... so its D=0
        TopD,                           // Lower Base  D
        DeflConeH,
        Density(0.0)                    // No Propellant there
      );

    // "EquipBay" is assumed to be a Cylinder, its UpperBase is the LowBase of
    // "DeflectorCone", ie @ TopX:
    constexpr static TC EquipBayProto =
      TC
      (
        TopX,
        TopD,
        EquipBayH,
        Density(0.0)                    // No propellant there
      );

    // The TrCone enclosing the OxidTankTop (which is a SpherSegm). Its Upper-
    // Base is LowBase of "EquipBay", so its D is "TopD" as well:
    //
    constexpr static TC  OxidTankTopEnclProto =
      TC
      (
        EquipBayProto.GetLow()[0],
        TopD,
        OxidTankTopD,
        OxidTankTopH,
        Density(0.0)                    // No propellant in it by itself...
      );

    //-----------------------------------------------------------------------//
    // OxidTank Protos:                                                      //
    //-----------------------------------------------------------------------//
    // The Top SpherSegm of the OxidTank. Its UpperBase is the LowBase of the
    // above "OxidTankTopEnclProto":
    //
    constexpr static SpS OxidTankTopProto =
      SpS
      (
        true,                           // Yes, Facing Up!
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

    constexpr static TC OxidTankUpProto =
      TC
      (
        OxidTankTopProto.GetLow()[0],
        OxidTankTopD,
        MaxD,
        OxidTankUpH,
        Propellants::LOxDens
      );

    // The Lower TrCone of the OxidTank. Its UpperBase is the LowerBase of the
    // "OxidTankUpProto":
    constexpr static TC OxidTankLowProto =
      TC
      (
        OxidTankUpProto.GetLow()[0],
        MaxD,
        OxidTankBtmD,
        OxidTankLowH,
        Propellants::LOxDens
      );

    // The Bottom SpherSegm of the OxidTank. Its Base is the LowBase of the
    // "OxidTankLowProto":
    constexpr static SpS OxidTankBtmProto =
      SpS
      (
        false,                      // Facing Low!
        OxidTankLowProto.GetLow()[0],
        OxidTankBtmD,
        OxidTankBtmH,
        Propellants::LOxDens
      );
/*
  public:
    //-----------------------------------------------------------------------//
    // "MechElement"s with Real Masses:                                      //
    //-----------------------------------------------------------------------//
*/

  public:
    //=======================================================================//
    // Thrust Vector Control:                                                //
    //=======================================================================//
    // Achived via 4 Gimbaled Vernier Chambers. Each Vernier Chamber is located
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
