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

  public:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Stage1 params are contained in the Base:
    //-----------------------------------------------------------------------//
    // Const LV and Mission Params:                                          //
    //-----------------------------------------------------------------------//
    // Number of Engines in Stage1:
    constexpr static int    NE = 9;

    // Engines Participating the the Boost-Back Burn:
    constexpr static double BBBurnEngPart    = 3.0 / double(NE);

    // Engines Participating in the Entry Burn:
    constexpr static double EntryBurnEngPart = 3.0 / double(NE);

    // Engines Participating in the Landing Burn:
    constexpr static double LandBurnEngPart  = 1.0 / double(NE);

    // NB: In this class, t=0 corresponds to the Stage1 Separation event, and
    // "t" runs FORWARD.
    // We assume that the landing site is (h=0, l=0). The co-ords at Separation
    // are (hS > 0, lS > 0), and they are considered to be const params:
    LenK   const          m_hS;
    LenK   const          m_lS;

    // And the corresp Velocity and Trajectory Inclidation at Separation:
    VelK   const          m_VS;
    Angle  const          m_psiS;

  private:
    // Params of a "Fully-Prop-Loaded" Stage1 (nominal):
    Mass   const          m_fplMass1;
    double const          m_fplK1;
    double const          m_fplPropRem1;
    Len    const          m_diam;

    //-----------------------------------------------------------------------//
    // Optimisation Params:                                                  //
    //-----------------------------------------------------------------------//
    // "PropMassS" is the propellant mass at Stage1 Separation, used to achieve
    // the soft RTLS and landing. Is to be MINIMISED subject to all constraints.
    // No need for a separate fld for it -- it will become Base::m_propMass1 !

    // Coast Dur (before BBBurn): Equals to BBIgnTime:
    Time                  m_coastDur;

    // Boost-Back Burn: Throttling Level decreases from "m_bbBurnThrtL0"  to
    // "m_bbBurnThrtL1" linearly in time,    and the Thrust vector elevation
    // changes linearly from "m_bbBurnTheta0" to "m_bbBurnTheta1"; both angles
    // should be around Pi (with cos(theta) < 0), obviously:
    Time                  m_bbBurnDur;
    double                m_bbBurnThrtL0;
    double                m_bbBurnThrtL1;
    Angle                 m_bbBurnTheta0;
    Angle                 m_bbBurnTheta1;

    // The trigger for the Entry (Slowing-Down) Burn is based on the Q, not the
    // temporal separation. Its purpose it to limit the Q durting re-entry.  We
    // assume that this burn occurs at full-thrust, with AoA=0, so the only var-
    // iable param is the burn duration.   We again assume a linear decrease of
    // the throttling level:
    Pressure              m_entryBurnQ;
    Time                  m_entryBurnDur;
    double                m_entryBurnThrtL0;
    double                m_entryBurnThrtL1;

    // XXX: For the Final Decsent and Landing, we currently do NOT perform any
    // special maneuvers  to avoid the "ballistic target" point and fly to the
    // actual landing site. We just perform a continuous burn in order to land
    // with near-zero velocity. So we currently assume AoA=0 and a const  (but
    // subject to optimisation) "mu", given by the corresp ThrottlingLevel.  A
    // reasonable trigger is the LandBurn altitude; the burn lasts until land-
    // ing (h=0), until v=0 (which means unsuccessful langing if "h" is above
    // the threshold), or until the propellant is exhausted.  We again assume
    // a linear decrease in the throttling level:
    LenK                  m_landBurnH;
    double                m_landBurnThrtL0;
    double                m_landBurnThrtL1;

    // So altogether: 14 params:
    // [
    //  PropMassS,  CoastDur,
    //  BBBurnDur,  BBBurnThrtL0,   BBBurnThrtL1,    BBBurnTheta0, BBBurnTheta1,
    //  EntryBurnQ, EntryBurnDur,   EntryBurnThrtL0, EntryBurnThrtL1,
    //  LandBurnH,  LandBurnThrtL0, LandBurnThrtL1
    // ]:
    constexpr static int  NP = 14;

    //-----------------------------------------------------------------------//
    // Scaling Factors for Translation of Rel Opt Params into the Abs Ones:  //
    //-----------------------------------------------------------------------//
    constexpr static Mass MaxPropMassS    = 50000.0_kg;
    constexpr static Time MaxBBBurnDur    =    50.0_sec;
    constexpr static Time MaxEntryBurnDur =    30.0_sec;
    constexpr static LenK MaxLandBurnH    =     3.0_km;

    //-----------------------------------------------------------------------//
    // Transient Data (during flight path integration):                      //
    //-----------------------------------------------------------------------//
    FlightMode            m_mode;
    Time                  m_entryIgnTime;   // Triggered by Q
    Time                  m_landIgnTime;    // Triggered by H
    Time                  m_maxLandBurnDur; // Auto-computed
    std::string           m_eventStr;
    Time                  m_nextOutputTime;

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
      Mass           a_fpl_mass1,           // Assuming Full-Prop-Load
      double         a_fpl_K1,              // Based on "a_fpl_mass1"
      double         a_fpl_prop_rem1,       // ditto
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
      Angle          a_psiS,

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
    // Ctor from a "Proto" and the Normalised Ctl Params:                    //
    //-----------------------------------------------------------------------//
    // So altogether: 14 params:
    // [
    //  PropMassS,  CoastDur,
    //  BBBurnDur,  BBBurnThrtL0,   BBBurnThrtL1,    BBBurnTheta0, BBBurnTheta1,
    //  EntryBurnQ, EntryBurnDur,   EntryBurnThrtL0, EntryBurnThrtL1,
    //  LandBurnH,  LandBurnThrtL0, LandBurnThrtL1
    // ]:

    RTLS1
    (
      RTLS1  const&               a_proto,
      Pressure                    a_Q_limit,
      double a_propMassSN,        double a_coastDurN,
      double a_bbBurnDurN,        double a_bbBurnThrtL0N,
      double a_bbBurnThrtL1N,     double a_bbBurnTheta0N,
      double a_bbBurnTheta1N,
      double a_entryBurnQN,       double a_entryBurnDurN,
      double a_entryBurnThrtL0N,  double a_entryBurnThrtL1N,
      double a_landBurnHN,        double a_landBurnThrtL0N,
      double a_landBurnThrtL1N
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
    bool ODECB(StateV* a_s, Time a_t, Time a_dt);

    // "AeroDynForces":
    // Atmospheric Conditions and Aerodynamic Drag Force:
    //
    std::tuple<EAM::AtmConds, ForceK, ForceK>
    AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const;

    // "PropBurnRate": (may be variable over time, >= 0):
    MassRate PropBurnRate(Mass a_m, Time a_t) const;

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
      // Data Flds: Same NP params:
      Mass      const m_propMassS;
      Time      const m_coastDur;
      Time      const m_bbBurnDur;
      double    const m_bbBurnThrtL0;
      double    const m_bbBurnThrtL1;
      Angle     const m_bbBurnTheta0;
      Angle     const m_bbBurnTheta1;
      Pressure  const m_entryBurnQ;
      Time      const m_entryBurnDur;
      double    const m_entryBurnThrtL0;
      double    const m_entryBurnThrtL1;
      LenK      const m_landBurnH;
      double    const m_landBurnThrtL0;
      double    const m_landBurnThrtL1;

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
      // LV Proto:
      RTLS1 const*                a_proto,
      // Optimisation Params:
      std::vector<double>*        a_init_vals,
      double                      a_min_prop_massSN,
      double                      a_min_bbb_durN,
      LenK                        a_land_dL_limit,
      VelK                        a_land_V_limit,
      double                      a_land_accG_limit, // Final Acc in "g" units
      Pressure                    a_Q_limit,
      double                      a_longG_limit,
      // NOMAD Params:
      int                         a_max_evals,
      int                         a_opt_seed,
      bool                        a_stop_if_feasible,
      double                      a_use_vns,         // In [0..1), 0: no VNS
      bool                        a_use_mt           // true unless debugging
    );
  };
}
// End namespace SpaceBallistics
