// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Missions/Ascent2.h":                  //
//                  Ascent-to-Orbit for a "Model" 2-Stage LV                 //
//===========================================================================//
#include "SpaceBallistics/Missions/LVBase.h"
#include <boost/property_tree/ptree.hpp>

namespace SpaceBallistics
{
  namespace EAM = EarthAtmosphereModel;

  //=========================================================================//
  // "Ascent2" Class:                                                        //
  //=========================================================================//
  class Ascent2: public LVBase<Ascent2>
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // For Constrained Optimisation:
    // The max number of params used in Optimisation: 13:
    // [thrustMult2, bHat2, muHat2, aAoAHat2, bAoAHat2, TGapRel,
    //  thrustMult1, bHat1, muHat1, aAoAHat1, bAoAHat1, alpha1, payLoadMassRel]:
    constexpr static int      NP             = 13;

    // Atmospheric Pressure at which Fairing Separation occurs (XXX: should it
    // be made configurable?):
    constexpr static Pressure FairingSepCond = Pressure(1.0);
    constexpr static Time     MaxTGap        = 1000.0_sec;      // Very Large!

    // Constraints on Start Altitude and Start Velocity (both should be close
    // to 0):
    constexpr static LenK     MaxStartH      = 0.1_km;          // 100 m
    constexpr static VelK     MaxStartV      = VelK(0.03);      //  30 m/sec

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The Base Class:                                                       //
    //-----------------------------------------------------------------------//
    using  Base = LVBase<Ascent2>;
    friend class  LVBase<Ascent2>;

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
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Over-All:                                                             //
    //-----------------------------------------------------------------------//
    // Const Over-All Params (not affected by Optimisation):
    Mass     const        m_maxStartMass;
    Mass     const        m_fairingMass;

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

    // Non-Const Stage2 Params (may be updated in the course of Optimisation):
    Mass                  m_fullMass2;
    Mass                  m_emptyMass2;
    Mass                  m_propMass2;
    Mass                  m_unSpendable2;
    Mass                  m_spendable2;
    ForceK                m_thrustVacI2; // Nominal (@ Max BurnRate)
    MassRate              m_burnRateI2;  // BurnRate @ Stage2 IgnTime
    Time                  m_T2;          // Nominal Stage2 BurnTime

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
    double                m_thrustMult2;
    MassT3                m_aMu2;
    MassT2                m_bMu2;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_bHat2;
    double                m_muHat2;

    // AoA for Stage2: AoA(t) = t * (a * t - b),
    // where t <= 0 is our standard integration backward-running time    (t=0
    // corresponds  to the orbital insertion instant, and AoA=0 at that point);
    // NB: the (-b) coeff !!!
    Angle    const        m_maxAoA2;
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
    double                m_thrustMult1;
    MassT3                m_aMu1;
    MassT2                m_bMu1;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_bHat1;
    double                m_muHat1;

    // AoA for Stage1: for similarity with Stage2,
    // AoA(tau) = tau * (a * tau + b),
    // where tau  = t - tIgn1 >= 0, so AoA=0 @ tau=0 (Stage1 ignition):
    Angle    const        m_maxAoA1;
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

      // Integration / Output Params:
      Time            a_ode_integr_step,
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
    bool ODECB(StateV* a_s, Time a_t, Time a_tau);

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
    class        NOMADEvaluator;
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

  private:
    //=======================================================================//
    // "RunNOMAD":                                                           //
    //=======================================================================//
    // Helper invoked from "FindOptimalAscentCtls". Returns "true" on success
    // (then "a_init_vals" contains the optimal params),   "false" otherwise:
    //
    static bool RunNOMAD
    (
      Ascent2*                            a_proto,
      bool                const           a_act_opts[NP],
      std::vector<double>*                a_init_vals,
      std::vector<double> const&          a_lo_bounds,
      std::vector<double> const&          a_up_bounds,
      int                                 a_max_evals,
      bool                                a_constr_q,
      bool                                a_constr_sep_q,
      bool                                a_constr_long_g,
      boost::property_tree::ptree const&  a_pt
    );
  };
}
// End namespace SpaceBallistics
