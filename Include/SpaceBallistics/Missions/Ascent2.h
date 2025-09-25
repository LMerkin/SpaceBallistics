// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Missions/Ascent2.h":                  //
//                  Ascent-to-Orbit for a "Model" 2-Stage LV                 //
//===========================================================================//
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include <boost/functional/hash.hpp>
#include <boost/property_tree/ptree.hpp>
#include <ostream>
#include <optional>
#include <utility>
#include <unordered_map>

// XXX: NOMAD headers produce tons of warnings, suppress them all:
#ifdef   __clang__     
#pragma  clang diagnostic push
#pragma  clang diagnostic ignored "-Wunused-parameter"
#pragma  clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma  clang diagnostic ignored "-Wcast-qual"
#pragma  clang diagnostic ignored "-Wcast-align"
#pragma  clang diagnostic ignored "-Wold-style-cast"
#pragma  clang diagnostic ignored "-Wsuggest-override"
#pragma  clang diagnostic ignored "-Wsuggest-destructor-override"
#pragma  clang diagnostic ignored "-Woverloaded-virtual" 
#pragma  clang diagnostic ignored "-Wshadow"
#pragma  clang diagnostic ignored "-Wextra-semi"
#pragma  clang diagnostic ignored "-Wunused-member-function"
#pragma  clang diagnostic ignored "-Winconsistent-missing-destructor-override"
#pragma  clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#pragma  clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma  clang diagnostic ignored "-Wdeprecated-dynamic-exception-spec"
#pragma  clang diagnostic ignored "-Wreserved-identifier"
#pragma  clang diagnostic ignored "-Wheader-hygiene"
#pragma  clang diagnostic ignored "-Wsign-conversion"
#pragma  clang diagnostic ignored "-Warray-bounds"
#pragma  clang diagnostic ignored "-Wmissing-noreturn"
#else
#pragma  GCC   diagnostic push
#pragma  GCC   diagnostic ignored "-Wunused-parameter"
#pragma  GCC   diagnostic ignored "-Woverloaded-virtual="
#pragma  GCC   diagnostic ignored "-Wnon-virtual-dtor"
#pragma  GCC   diagnostic ignored "-Wcast-qual"
#pragma  GCC   diagnostic ignored "-Warray-bounds"
#endif
#include <Nomad/nomad.hpp>
#include <Cache/CacheBase.hpp>
#ifdef   __clang__
#pragma  clang diagnostic pop
#else
#pragma  GCC   diagnostic pop
#endif

namespace SpaceBallistics
{
  namespace EAM = EarthAtmosphereModel;

  //=========================================================================//
  // The "Ascent2" Class:                                                    //
  //=========================================================================//
  class Ascent2
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    constexpr static GMK      K              = BodyData<Body::Earth>::K;
    constexpr static LenK     R              = 6371.0_km;   // Equi-Volume

    // ODE Integration Params: 1 msec step; it may only be reduced, never
    // increased beyond the original value:
    constexpr static Time     ODEInitStep    = 0.001_sec;
    constexpr static Time     ODEMaxStep     = ODEInitStep;
    constexpr static double   ODERelPrec     = 1e-6;

    // Singular point detection criteria: NB: "Vhor" approaches 0 much faster
    // than "Vr":
    constexpr static VelK     SingVr         = VelK(1e-2); // 10   m/sec
    constexpr static VelK     SingVhor       = VelK(1e-4); //  0.1 m/sec

    // For Constrained Optimisation:
    // The max number of params used in Optimisation: 13:
    // [thrustMult2, bHat2, muHat2, aAoAHat2, bAoAHat2, TGapRel,
    //  thrustMult1, bHat1, muHat1, aAoAHat1, bAoAHat1, alpha1, payLoadMassRel]:
    constexpr static int      NP             = 13;

    // Atmospheric Pressure at which Fairing Separation occurs (XXX: should it
    // be made configurable?):
    constexpr static Pressure FairingSepCond = Pressure(1.0);

    // Constraints on Start Altitude and Start Velocity (both should be close
    // to 0):
    constexpr static Time     MaxTGap        = 1000.0_sec;      // Very Large!
    constexpr static LenK     MaxStartH      = 0.1_km;          // 100 m
    constexpr static VelK     MaxStartV      = VelK(0.03);      //  30 m/sec

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    using MassT2 =   decltype(MassRate(1.0) / 1.0_sec);
    using MassT3 =   decltype(MassT2  (1.0) / 1.0_sec);

    //-----------------------------------------------------------------------//
    // "StateV":                                                             //
    //-----------------------------------------------------------------------//
    // State Vector (for some time t <= 0, where t=t0=0 corresponds to the Orbi-
    // tal Insertion):
    // Computations are performed in TopoCentric (Start) Planar Polar CoOrds, so
    // the first 3 components are r, rDot, omega = phiDot. The 4th  component is
    // the Total Spent Propellant Mass between the curr "t" and "t0" (which is a
    // continuous and DECREASING function of "t",  as opposed to the Total Mass
    // which is discontinuous when the Stage1 or Fairing are jettisones), and
    // the 5th component is the polar angle "phi" (integrated omega):
    //
    using StateV  = std::tuple<LenK, VelK, AngVel, Mass, Angle>;
    //                         r   rDot=Vr omega   spent phi

    //-----------------------------------------------------------------------//
    // "DStateV":                                                            //
    //-----------------------------------------------------------------------//
    // The Time Derivative of the "StateV". The "MassRate" components is the
    // negated (<= 0) BurnRate:
    //
    using DStateV = std::tuple<VelK, AccK, AngAcc, MassRate, AngVel>;

    //-----------------------------------------------------------------------//
    // "FlightMode":                                                         //
    //-----------------------------------------------------------------------//
    enum class FlightMode: int
    {
      UNDEFINED = 0,
      Burn1     = 1,  // Stage1 Burn,
      Gap       = 2,  // Ballistic Gap
      Burn2     = 3   // Stage2 Burn
    };

    static char const* ToString(FlightMode a_mode)
    {
      switch (a_mode)
      {
        case FlightMode::Burn1: return "Burn1";
        case FlightMode::Gap  : return "Gap";
        case FlightMode::Burn2: return "Burn2";
        default:                return "UNDEFINED";
      }
    }

  private:
    //-----------------------------------------------------------------------//
    // "NearSingularityException" Class:                                     //
    //-----------------------------------------------------------------------//
    // Exception thrown when the flight path is approaching the singular point:
    //
    class NearSingularityExn
    {
    public:
      // Data Flds:
      LenK  const m_r;
      VelK  const m_Vr;
      Mass  const m_spentPropMass;
      Angle const m_phi;
      Time  const m_t;

      // Non-Default Ctor:
      NearSingularityExn(LenK  a_r,   VelK a_vr, Mass a_spent_prop_mass,
                         Angle a_phi, Time a_t)
      : m_r             (a_r),
        m_Vr            (a_vr),
        m_spentPropMass (a_spent_prop_mass),
        m_phi           (a_phi),
        m_t             (a_t)
      {}
    };

  public:
    //-----------------------------------------------------------------------//
    // "RunRC":                                                              //
    //-----------------------------------------------------------------------//
    // Possible results of "Run":
    //
    enum class RunRC: int
    {
      // The "desirable" result: Singular Point reached: V = 0 at some h >= 0:
      Singularity = 0,

      // Reached h = 0 at V > 0, still with unspent propellant (beyond the
      // minimum remnant):
      ZeroH       = 1,

      // Ran out of all available propellant (up to the minimum remnant), but
      // still h > 0 and V > 0:
      FlameOut    = 2,

      // If any exception occurred:
      Error       = 3
    };

    static char const* ToString(RunRC a_rc)
    {
      switch (a_rc)
      {
        case RunRC::Singularity: return "Singularity";
        case RunRC::ZeroH      : return "ZeroH";
        case RunRC::FlameOut   : return "FlameOut";
        default                : return "Error";
      }
    }

    //-----------------------------------------------------------------------//
    // "RunRes" Struct:                                                      //
    //-----------------------------------------------------------------------//
    struct RunRes
    {
      RunRC     m_rc;       // Return Code
      Time      m_T;        // Final (actually start) Time, < 0
      LenK      m_hT;       // Final (actually start) Altitude
      VelK      m_VT;       // Final (sctually start) Velocity
      Mass      m_mT;       // Final (actually start) Mass
      Pressure  m_maxQ;     // Max Dynamic Pressure (Q)
      Pressure  m_sepQ;     // Q @ Stage1 Separation
      double    m_maxLongG; // Max Longitudinal G
    };

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Over-All:                                                             //
    //-----------------------------------------------------------------------//
    // Const Over-All Params (not affected by Optimisation):
    using    AreaK = decltype(Sqr(1.0_km));
    Mass     const        m_maxStartMass;
    Mass     const        m_fairingMass;
    AreaK    const        m_crosS;

    // Non-Const Over-All Params (may be updated in the course of Optimisation):
    double                m_alpha1;
    Mass                  m_payLoadMass;

    //-----------------------------------------------------------------------//
    // Stage2:                                                               //
    //-----------------------------------------------------------------------//
    // Const Stage2 Params (not affected by Optimisation):
    double   const        m_K2;
    double   const        m_propRem2;
    Time     const        m_IspVac2;
    double   const        m_minThrtL2;
    Angle    const        m_maxAoA2;

    // Non-Const Stage2 Params (may be updated in the course of Optimisation):
    Mass                  m_fullMass2;
    Mass                  m_emptyMass2;
    Mass                  m_propMass2;
    Mass                  m_unSpendable2;
    Mass                  m_spendable2;
    ForceK                m_thrustVacI2; // Nominal (@ Max BurnRate)
    double                m_thrustMult2;
    MassRate              m_burnRateI2;  // BurnRate @ Stage2 IgnTime
    Time                  m_T2;          // Nominal Stage2 BurnTime

    //-----------------------------------------------------------------------//
    // Stage1:                                                               //
    //-----------------------------------------------------------------------//
    // Const Stage1 Params (not affected by Optimisation):
    double   const        m_K1;
    double   const        m_propRem1;
    Time     const        m_IspSL1;
    Time     const        m_IspVac1;
    double   const        m_minThrtL1;
    Angle    const        m_maxAoA1;

    // Non-Const Stage1 Params (may be updated in the course of Optimisation):
    Mass                  m_fullMass1;
    Mass                  m_emptyMass1;
    Mass                  m_propMass1;
    Mass                  m_unSpendable1;
    Mass                  m_spendable1;
    ForceK                m_thrustVacI1; // Nominal (@ Max BurnRate)
    double                m_thrustMult1;
    MassRate              m_burnRateI1;  // BurnRate @ Stage1 IgnTime
    Time                  m_T1;          // Nominal Stage1 BurnTime

    //-----------------------------------------------------------------------//
    // Mission Params:                                                       //
    //-----------------------------------------------------------------------//
    LenK                  m_Rins;        // Radius-Vector @ Orbital Insertion
    VelK                  m_Vins;        // LV Velocity   @ Orbital Insertion

    //-----------------------------------------------------------------------//
    // Flight Control Program Parameterisation:                              //
    //-----------------------------------------------------------------------//
    // BurnRate for Stage2:
    // It is a non-increasing quadratic function of time:
    // BurnRate(tau) = BurnRateI + bMu * tau + aMu * tau^2, where "tau" is the
    // time since Stage2 ignition (0 <= tau <= T2):
    MassT3                m_aMu2;
    MassT2                m_bMu2;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_bHat2;
    double                m_muHat2;

    // AoA for Stage2: AoA(t) = t * (a * t - b),
    // where t <= 0 is our standard integration backward-running time    (t=0
    // corresponds  to the orbital insertion instant, and AoA=0 at that point);
    // NB: the (-b) coeff !!!
    AngAcc                m_aAoA2;
    AngVel                m_bAoA2;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_aAoAHat2;
    double                m_bAoAHat2;

    // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
    // ignition):
    Time                  m_TGap;

    // BurnRate for Stage1: Similar to that of Stage1, where "tau" is the time
    // since Stage1 ignition (0 <= tau <= T1):
    MassT3                m_aMu1;
    MassT2                m_bMu1;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_bHat1;
    double                m_muHat1;

    // AoA for Stage1: for similarity with Stage2,
    // AoA(tau) = tau * (a * tau + b),
    // where tau  = t - tIgn1 >= 0, so AoA=0 @ tau=0 (Stage1 ignition):
    AngAcc                m_aAoA1;
    AngVel                m_bAoA1;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_aAoAHat1;
    double                m_bAoAHat1;

    //-----------------------------------------------------------------------//
    // Constraints used in Optimisation:                                     //
    //-----------------------------------------------------------------------//
    // (see below for the actual encountered vals):
    Pressure              m_QLimit;
    Pressure              m_sepQLimit;
    double                m_longGLimit;

    //-----------------------------------------------------------------------//
    // Transient Data (during flight path integration):                      //
    //-----------------------------------------------------------------------//
    FlightMode            m_mode;
    Time                  m_ignTime2;       // Burn2 start: St2 Ignition Time
    Time                  m_fairingSepTime; // Fairing Separation Time
    Time                  m_cutOffTime1;    // Gap   start: St1 Cut-Off  Time
    Time                  m_ignTime1;       // Burn1 start: St1 Ignition Time

    // Vals to be Constrined:
    Pressure              m_maxQ;           // Max Dynamic Pressure (Q)
    Pressure              m_sepQ;           // Q @ Stage1 separation
    double                m_maxLongG;       // Max Longitudinal G

    // For output:
    std::ostream*         m_os;
    int                   m_logLevel;

  public:
    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    Ascent2
    (
      // Stage2:
      double          a_K2,             // PropMass2 / FullMass2
      double          a_prop_rem2,
      Time            a_Isp_vac2,
      ForceK          a_thrust_vac2,
      double          a_min_thrtl2,
      Angle_deg       a_max_aoa2,

      // Stage1:
      double          a_K1,             // PropMass1 / FullMass1
      double          a_prop_rem1,
      Time            a_Isp_sl1,
      Time            a_Isp_vac1,
      ForceK          a_thrust_vac1,
      double          a_min_thrtl1,
      Angle_deg       a_max_aoa1,

      // Over-All:
      double          a_alpha1,         // FullMass1 / FullMass2
      Mass            a_max_start_mass,
      Mass            a_fairing_mass,
      Len             a_diam,
      Mass            a_payload_mass,

      // Constraints (+oo if no constraint):
      Pressure        a_Q_limit,
      Pressure        a_sepQ_limit,
      double          a_longG_limit,

      // Mission Params:
      LenK            a_h_perigee,
      LenK            a_h_apogee,
      Angle_deg       a_incl,
      Angle_deg       a_launch_lat,

      // Logging Params:
      std::ostream*   a_os,
      int             a_log_level
    );

    //-----------------------------------------------------------------------//
    // Copy Ctor: Required for optimisation:                                 //
    //-----------------------------------------------------------------------//
    Ascent2(Ascent2 const&) = default;

    //-----------------------------------------------------------------------//
    // "SetCtlParams":                                                       //
    //-----------------------------------------------------------------------//
    // Sets the Ctl Params which are not set by the above Ctor (to avoid over-
    // complexity of the latter):
    //
    void SetCtlParams
    (
      double   a_bHat2,    double   a_muHat2,
      double   a_aAoAHat2, double   a_bAoAHat2,
      Time     a_TGap,
      double   a_bHat1,    double   a_muHat1,
      double   a_aAoAHat1, double   a_bAoAHat1
    );

    //-----------------------------------------------------------------------//
    // "Run": Integrate the Ascent Trajectory:                               //
    //-----------------------------------------------------------------------//
    RunRes Run();

  private:
    //=======================================================================//
    // Internal Helpers:                                                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // For ODE Integration:                                                  //
    //-----------------------------------------------------------------------//
    // "ODERHS":
    // For Ascent Trajectory Integration:
    //
    DStateV ODERHS(StateV const& a_s, Time a_t);

    // "ODECB":
    // CallBack (invoked after the completion of each integration step):
    // Here FlightMode switching occurs, so this method is non-"const":
    //
    bool ODECB(StateV* a_s, Time a_t);

    // Common part of "ODERHS" and "ODECB":
    void NonGravForces
    ( 
      // Inputs:
      Time           a_t,
      LenK           a_r,
      VelK           a_Vr,
      VelK           a_Vhor,
      Mass           a_m,
      // Outputs:
      VelK*          a_V,
      Angle*         a_psi,
      Angle*         a_aoa,
      EAM::AtmConds* a_atm, 
      MassRate*      a_burn_rate,
      ForceK*        a_thrust,
      double         a_lv_axis[2],  // In the (r, normal-to-r) frame
      AccK           a_ng_acc [2]   // ditto
    ) 
    const;

    // "AeroDynForces":
    // Atmospheric Conditions and Aerodynamic Drag Force:
    //
    std::tuple<EAM::AtmConds, ForceK, ForceK>
    AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const;

    // "LocateSingularPoint":
    // The final stage of integration, where we assume omega=0 (a purely vert-
    // ical motion)   and integrate the simplified ODEs ANALYTICALLY to avoid
    // numerical instabilities in the viciniy of the singular point.
    // Returns (singH, singT)  if the singular point has been found, otherwise
    // "nullopt":
    std::optional<std::pair<StateV, Time>>
    LocateSingularPoint(NearSingularityExn const& a_nse);

    //      "PropBurnRate": (may be variable over time, >= 0):
    MassRate PropBurnRate(Time a_t) const;

    //   "AoA": Angle-of-Attack (variable over time, also constrained with the
    //   curr pitch "psi"):
    Angle AoA(Time a_t, Angle a_psi) const;

    //    "Thrust": Depends on the Mode, BurnRate and the Counter-Pressure:
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const;

    //  "LVMass": Current Mass (LV + PayLoad):
    Mass LVMass(StateV const& a_s, Time a_t) const;

    //  "OutputCtls": For Testing Only:
    void OutputCtls() const;

    //-----------------------------------------------------------------------//
    // For Optimisation:                                                     //
    //-----------------------------------------------------------------------//
    // NB: "ModifyLVParams" also includes "SetCtlParams":
    //
    void ModifyLVParams
    (
      // Thrust and Mass Params:
      double a_thrustMult2, double a_thrustMult1,
      double a_alpha1,      Mass   a_payLoadMass,
      // Ctl Params:
      double a_bHat2,       double a_muHat2,
      double a_aAoAHat2,    double a_bAoAHat2,
      Time   a_TGap,
      double a_bHat1,       double a_muHat1,
      double a_aAoAHat1,    double a_bAoAHat1
    );

    //=======================================================================//
    // "NOMADEvaluator": Helper Class used in NOMAD Optimisation:            //
    //=======================================================================//
    class NOMADEvaluator final: public NOMAD::Evaluator
    {
    private:
      //---------------------------------------------------------------------//
      // Data Flds:                                                          //
      //---------------------------------------------------------------------//
      Ascent2 const*  m_proto;
      bool            m_actOpts[NP];    // Flags for Active Optimisation Params

    public:
      //---------------------------------------------------------------------//
      // Non-Default Ctor, Dtor:                                             //
      //---------------------------------------------------------------------//
      NOMADEvaluator
      (
        Ascent2 const*                                a_proto,
        bool    const                                 a_actOpts[NP],
        std::shared_ptr<NOMAD::EvalParameters> const& a_params
      );

      ~NOMADEvaluator() override {}

      //---------------------------------------------------------------------//
      // "eval_x": The actual evaluation method for NOMAD:                   //
      //---------------------------------------------------------------------//
      bool eval_x
      (
        NOMAD::EvalPoint&    a_x,  // Not "const" -- the result is also set here
        NOMAD::Double const& UNUSED_PARAM(a_hMax),
        bool&                a_countEval
      )
      const override;
    };
    friend class NOMADEvaluator;

  public:
    //=======================================================================//
    // "OptRes" Struct:                                                      //
    //=======================================================================//
    // The result of "FindOptimalAscentCtls" below. It just contains the rele-
    // vant flds taken from the main "Ascent2" class:
    //
    struct OptRes
    {
      //=====================================================================//
      // Data Flds:                                                          //
      //=====================================================================//
      //---------------------------------------------------------------------//
      // Stage2 Params:                                                      //
      //---------------------------------------------------------------------//
      // Stage2 Thrust Multiplier and Thrust (@ Ignition):
      double    m_thrustMult2;
      ForceK    m_thrustVacI2;
      MassRate  m_burnRateI2;

      // Stage2 BurnRate: Dim-Less Params:
      double    m_bHat2;            // Must be in [0 .. 1], default is 0
      double    m_muHat2;           // Must be in [0 .. 1], default is 1

      // Stage2 BurnRate: Dimensioned Coeffs:
      Time      m_T2;
      MassT3    m_aMu2;
      MassT2    m_bMu2;

      // Stage2 AoA: Dim-Less Params:
      double    m_aAoAHat2;         // Must be in [0 .. 1], default is 0
      double    m_bAoAHat2;         // Must be in [0 .. 1], default is 0

      // Stage2 AoA: Dimensioned Coeffs: AoA2(t) = t * (a2 * t - b2), t <= 0,
      // t=0 is the orbital insertion time:
      AngAcc    m_aAoA2;
      AngVel    m_bAoA2;

      //---------------------------------------------------------------------//
      // Ballistic Gap:                                                      //
      //---------------------------------------------------------------------//
      Time      m_TGap;             // Default is 0.0_sec

      //---------------------------------------------------------------------//
      // Stage1 Params:                                                      //
      //---------------------------------------------------------------------//
      // Stage1 Thrust Multiplier and Thrust (@ Ignition):
      double    m_thrustMult1;
      ForceK    m_thrustVacI1;
      MassRate  m_burnRateI1;

      // Stage1 BurnRate: Dim-Less Params:
      double    m_bHat1;            // Must be in [0 .. 1], default is 0
      double    m_muHat1;           // Must be in [0 .. 1], default is 1

      // Stage1 BurnRate: Dimensioned Coeffs:
      Time      m_T1;
      MassT3    m_aMu1;
      MassT2    m_bMu1;

      // Stage1 AoA: Dim-Less Params:
      double    m_aAoAHat1;         // Must be in [0 .. 1], default is 0
      double    m_bAoAHat1;         // Must be in [0 .. 1], default is 0

      // Stage1 AoA: Dimensioned Coeffs: AoA1(t) = tau * (a1 * tau + b1),
      // tau >= 0 is the time since Stage1 ignition:
      AngAcc    m_aAoA1;
      AngVel    m_bAoA1;

      //---------------------------------------------------------------------//
      // Over-All:                                                           //
      //---------------------------------------------------------------------//
      double    m_alpha1;           // FullMass1 / FullMass2
      Mass      m_payLoadMass;

      //=====================================================================//
      // Non-Default Ctor:                                                   //
      //=====================================================================//
      OptRes(Ascent2 const& a_asc);

      //=====================================================================//
      // Output:                                                             //
      //=====================================================================//
      friend std::ostream& operator<<
            (std::ostream& a_os, OptRes const& a_res);
    };
    friend struct OptRes;

    //=======================================================================//
    // "FindOptimalAscentCtls":                                              //
    //=======================================================================//
    // Top-Level Optiomisation Function.
    // Tries to find the Flight Control Params such that  the specified target
    // orbit is reached (with the orbital insertion occurring in the perigee,
    // if the orbit is elliptical) and that requires the min LV StartMass  (if
    // the PayloadMass is fixed)  or otherwise, achieves the max PayLoadMass.
    // Returns "OptRes" if successful, "nullopt" otherwise.
    // Can also perform the "final run" on the optimal params found, and return
    // the corresp "RunRes" (if specified in the Config.ini):
    //
    static std::pair<std::optional<OptRes>,
                     std::optional<RunRes>>
    FindOptimalAscentCtls
    (
      // All are are given via the ConfigFile.ini, as there are quite a few of
      // them:
      std::string const& a_config_ini,
      std::ostream*      a_os   // May be NULL
    );
  };
}
// End namespace SpaceBallistics
