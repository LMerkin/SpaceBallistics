// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Missions/RTLS1.h":                     //
//           Return-To-(Launch/Landing)-Site for Stage1 of an LV						 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Missions/LVBase.h"
#include <array>
#include <vector>
#include <utility>

namespace SpaceBallistics
{
  //=========================================================================//
  // "RTLS1" Class:                                                          //
  //=========================================================================//
  class RTLS1: public LVBase<RTLS1>
  {
  public:
    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    using  Base = LVBase<RTLS1>;
    friend class  LVBase<RTLS1>;

    // Number of Coeffs in sin(theta) expansion (so the degree is NS-1):
    constexpr static int  NS = 4;

  private:
    //-----------------------------------------------------------------------//
    // "FlightMode":                                                         //
    //-----------------------------------------------------------------------//
    enum class FlightMode: int
    {
      UNDEFINED   = 0,
      Coast       = 1, // From separation to BoostBackBurn (BBBurn): passive
      BBBurn      = 2, //
      ExoAtmDesc  = 3, // Exo-Atmospheric  Descent:                  passive
      EntryBurn   = 4, // Slowing-Down before Re-Entering the Atmosphere
      EndoAtmDesc = 5, // Endo-Atmospheric Descent:                  passive
      LandBurn    = 6  //
    };

    static char const* ToString(FlightMode a_mode)
    {
      switch (a_mode)
      {
        case FlightMode::Coast      : return "Coast";
        case FlightMode::BBBurn     : return "BBBurn";
        case FlightMode::ExoAtmDesc : return "ExoAtmDesc";
        case FlightMode::EntryBurn  : return "EntryBurn";
        case FlightMode::EndoAtmDesc: return "EndoAtmDesc";
        case FlightMode::LandBurn   : return "LandBurn";
        default:                      return "UNDEFINED";
      }
    }

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Stage1 params are contained in the Base:
    //-----------------------------------------------------------------------//
    // Const Mission Params:                                                 //
    //-----------------------------------------------------------------------//
    // NB: In this class, t=0 corresponds to the Stage1 Separation event, and
    // "t" runs FORWARD.
    // We assume that the landing site is (h=0, l=0). The co-ords at Separation
    // are (hS > 0, lS > 0), and they are considered to be const params:
    LenK  const           m_hS;
    LenK  const           m_lS;

    // And the corresp Velocity and Trajectory Inclidation at Separation:
    VelK  const           m_VS;
    Angle const           m_phiS;

    //-----------------------------------------------------------------------//
    // Optimisation Params:                                                  //
    //-----------------------------------------------------------------------//
    // "PropMassS" is the propellant mass at Stage1 Separation, used to achieve
    // the soft RTLS and landing. Is to be MINIMISED subject to all constraints.
    // No need for a separate fld for it -- it will become Base::m_propMass1 !

    // Coast Dur (before BBBurn): Equals to BBIgnTime:
    Time                  m_coastDur;

    // AoA ctl during BBBurn  (more precisely, we control "theta" -- the Thrust
    // vector inclination angle).   Generally, the Thrust vector points towards
    // (-x),  so "theta" is around Pi. We expand sin(theta) into a CUBIC (or lo-
    // wer degree) polynomial of (t_burn / burn_duration), and it can be of any
    // sign, whereas we keep cos(theta) < 0. There are no other restrictions on
    // "theta" / AoA in this case, because the BoostBackBurn is Exo-Atmospheric:    // XXX: Do we always need FullThrust here?
    Time                  m_bbBurnDur;
    double                m_bbBurnSinTheta[NS];

    // The trigger for the Entry (Slowing-Down) Burn is based on the Q, not the
    // temporal separation. Its purpose it to limit the Q durting re-entry.  We
    // assume that this burn occurs at full-thrust, with AoA=0, so the only var-
    // iable param is the burn duration:
    // XXX: Do we always need FullThrust here?
    Pressure              m_entryBurnQ;
    Time                  m_entryBurnDur;

    // XXX: For the Final Decsent and Landing, we currently do NOT perform any
    // special maneuvers  to avoid the "ballistic target" point and fly to the
    // actual landing site. We just perform a continuous burn in order to land
    // with near-zero velocity. So we currently assume AoA=0 and a const  (but
    // subject to optimisation) "mu", given by the corresp ThrottlingLevel.  A
    // reasonable trigger is the LandBurn altitude; the burn lasts until land-
    // ing (h=0), until v=0 (which means unsuccessful langing if "h" is above
    // the threshold), or until the propellant is exhausted:
    LenK                  m_landBurnH;
    double                m_landBurnThrtL;

    // So altogether: 11 params:
    // [PropMassS,  CoastTime,    BBBurnDur, BBBurnSinTheta[4],
    //  EntryBurnQ, EntryBurnDur, LandBurnH, LandBurnThrtL]
    constexpr static int  NP = 11;

    //-----------------------------------------------------------------------//
    // Translation of Relative Opt Params in the Absolute Ones:              //
    //-----------------------------------------------------------------------//
    // Scaling Factor:
    constexpr static Mass     MaxPropMassS     = 30000.0_kg;
    constexpr static Time     MaxCoastDur      =   300.0_sec;
    constexpr static Time     MaxBBBurnDur     =    20.0_sec;  // Too large?
    constexpr static Pressure MaxEntryBurnQ    = Pressure(30000.0);
    constexpr static Time     MaxEntryBurnDur  =    10.0_sec;  // Too large?
    constexpr static LenK     MaxLandBurnH     =     5.0_km;   // Too low?

    //-----------------------------------------------------------------------//
    // Transient Data (during flight path integration):                      //
    //-----------------------------------------------------------------------//
    FlightMode            m_mode;
    Time                  m_entryIgnTime;   // Triggered by Q
    Time                  m_landIgnTime;    // Triggered by H
    Time                  m_finalTime;
    std::string           m_eventStr;

    // Vals which may be Constrained:
    Pressure              m_maxQ;
    // The following is primarily for compatibility with the "LVBase":
    Pressure              m_sepQ;
    double                m_maxLongG;

    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    RTLS1
    (
      // Stage Params:
      Mass           a_full_mass1,          // Assuming full prop load
      double         a_full_K1,             // Based on "a_full_mass1"
      double         a_full_prop_rem1,      // ditto
      Time           a_Isp_sl1,
      Time           a_Isp_vac1,
      ForceK         a_thrust_vac1,
      double         a_min_thrtl1,
      Len            a_diam,

      // Mission (Return-To-Launch/Landing-Site) Params:
      Mass           a_prop_massS,          // INCL the unspendable remnant
      LenK           a_hS,
      LenK           a_lS,
      VelK           a_VS,
      Angle          a_phiS,

      // Integration and Output Params:
      Time           a_ode_integr_step,
      std::ostream*  a_os,
      int            a_log_level
    );

    //-----------------------------------------------------------------------//
    // Copy Ctor: Required for Optimisation:                                 //
    //-----------------------------------------------------------------------//
    RTLS1(RTLS1 const&) = default;

    //-----------------------------------------------------------------------//
    // "SetCtlParams":                                                       //
    //-----------------------------------------------------------------------//
    // See the implementation for details. There are 2 overloaded forms of this
    // methods.
    // IMPORTANT: "PropMassS" is NOT set by these methods; it is NOT a Ctl Par-
    // am; it can only be set when a new "RTLS1" obj is constructed!
    //
    // Vector Form:
    void SetCtlParams(std::vector<double> const& a_opt_params_n);

    // Individual Args Form:
    void SetCtlParams
    (
      double a_coastDurN,     double a_bbBurnDurN,
      double a_entryBurnQN,   double a_entryBurnDurN,
      double a_landBurnHN,    double a_landBurnThrtN,
      double a_sinTheta0,     double a_sinTheta1,
      double a_sinTheta2,     double a_sinTheta3
    );

    //-----------------------------------------------------------------------//
    // "Run": Integrate the Return Trajectory:                               //
    //-----------------------------------------------------------------------//
    Base::RunRes Run();

  private:
    //=======================================================================//
    // Internal Helpers:                                                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // For ODE Integration:                                                  //
    //-----------------------------------------------------------------------//
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

    // "PropBurnRate": (may be variable over time, >= 0):
    MassRate PropBurnRate(Time a_t) const;

    // "AoA": Angle-of-Attack (variable over time, also constrained with the
    // curr pitch "psi"):
    std::pair<Angle,Angle> AoA(Time a_t, Angle a_psi) const;

    // "Thrust": Depends on the Mode, BurnRate and the Counter-Pressure:
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const;

    // "LVMass": Current Mass:
    Mass LVMass(StateV const& a_s, Time a_t) const;

    //-----------------------------------------------------------------------//
    // "NOMADEvaluator": Helper Class used in NOMAD Optimisation:            //
    //-----------------------------------------------------------------------//
    class        NOMADEvaluator;
    friend class NOMADEvaluator;

  public:
    //-----------------------------------------------------------------------//
    // "OptRes" Struct: The result of "FindOptimalReturnCtls" below:         //
    //-----------------------------------------------------------------------//
    // It just contains the relevant flds taken from the main "RTLS1" class:
    //
    struct OptRes
    {
      // Data Flds:
      Mass      const m_propMassS;
      Time      const m_coastDur;
      Time      const m_bbBurnDur;
      double    const m_bbBurnSinTheta[NS];
      Pressure  const m_entryBurnQ;
      Time      const m_entryBurnDur;
      LenK      const m_landBurnH;
      double    const m_landBurnThrtL;

      // Non-Default Ctor:
      OptRes(RTLS1 const& a_rtls);

      // Output:
      friend std::ostream& operator<<
            (std::ostream& a_os, OptRes const& a_res);
    };
    friend struct OptRes;

    //-----------------------------------------------------------------------//
    // "FindOptimalReturnCtls":                                              //
    //-----------------------------------------------------------------------//
    // Top-Level Optimisation Function. Tries to find the ctls which minimise
    // the Propellant Mass and ensures all Return and Landing constraints are
    // satisfied. Can also perform the "final run" on the optimal params found,
    // and return the corresp "RunRes" (if specified in the Config.ini):
    //
    static std::pair<std::optional<OptRes>,
                     std::optional<RunRes>>
    FindOptimalReturnCtls
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
    // Helper invoked from "FindOptimalReturnCtls". Returns "true" on success
    // (then "a_init_vals" contains the optimal params),   "false" otherwise:
    //
    static bool RunNOMAD
    (
      // LV Params:
      RTLS1 const*                a_proto,
      Mass                        a_full_mass1,
      double                      a_full_k1,
      double                      a_full_prop_rem1,
      Len                         a_diam,
      // InitVals for the Ctl Params:
      std::vector<double>*        a_init_vals,
      // Optimisation Constraints (Limits):
      LenK                        a_land_dL_limit,
      VelK                        a_land_V_limit,
      Pressure                    a_Q_limit,
      // NOMAD Params:
      int                         a_max_evals,
      int                         a_opt_seed,
      bool                        a_stop_if_feasible,
      double                      a_use_vns     // In [0..1), 0: no VNS
    );
  };
}
// End namespace SpaceBallistics
