// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Missions/RTLS1.h":                     //
//                Return-To-Launch-Site for Stage1 of an LV									 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Missions/LVBase.h"
#include <array>

namespace SpaceBallistics
{
  //=========================================================================//
  // "RTLS1" Class:                                                          //
  //=========================================================================//
  class RTLS1: public LVBase<RTLS1>
  {
  private:
    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    using  Base = LVBase<Ascent2>;
    friend class  LVBase<Ascent2>;

    //-----------------------------------------------------------------------//
    // "FlightMode":                                                         //
    //-----------------------------------------------------------------------//
    enum class FlightMode: int
    {
      UNDEFINED = 0,
      Coast     = 1,  // From separation to BoostBackBurn (BBBurn): passive
      BBBurn    = 2,  //
      ExAtmDesc = 3,  // Exo-Atmospheric  Descent:                  passive
      EntryBurn = 4,  // Slowing-Down before Re-Entering the Atmosphere
      EnAtmDesc = 5,  // Endo-Atmospheric Descent:                  passive
      FinalBurn = 6   // 
    };

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    // Stage1 params are contained in the Base:
    using  Base = LVBase<RTLS1>;
    friend class  LVBase<RTLS1>;

    //-----------------------------------------------------------------------//
    // Const Mission Params:                                                 //
    //-----------------------------------------------------------------------//
    // NB: In this class, t=0 corresponds to the Stage1 Separation event, and
    // "t" runs FORWARD.
    // We assume that the landing site is (h=0, l=0). The co-ords at Separation
    // are (hS > 0, lS > 0), and they are considered to be const params:
    LenL const            m_hS;
    LenK const            m_lS;

    // And the corresp velocity components at Separation: Velocity  projections
    // on the radius-vector and on the positive normal to the radius-vector (ie
    // towards the increasing polar angle) are "Vr" and "Vhor", resp:
    VelK const            m_VrS;
    VelK const            m_VhorS;

    //-----------------------------------------------------------------------//
    // Optimisation Params:                                                  //
    //-----------------------------------------------------------------------//
    // "PropMassS" is the propellant mass at Stage1 Separation, used to achieve
    // the soft RTLS and landing. Is to be MINIMISED subject to all constraints:
    Mass                  m_propMassS;

    // Coast Time (before BBBurn):
    Time                  m_coastTime;

    // AoA ctl during BBBurn  (more precisely, we control "theta" -- the Thrust
    // vector inclination angle).   Generally, the Thrust vector points towards
    // (-x),  so "theta" is around Pi. We expand sin(theta) into a CUBIC (or lo-
    // wer degree) polynomial of (t_burn / burn_duration), and it can be of any
    // sign, whereas we keep cos(theta) < 0. There are no other restrictions on
    // "theta" / AoA in this case, because the BoostBackBurn is Exo-Atmospheric:    // XXX: Do we always need FullThrust here?
    Time                  m_bbBurnDur;
    constexpr static int  NS = 4;
    double                m_bbBurnSinTheta[NS];

    // The trigger for the Entry (Slowing-Down) Burn is based on the Q, not the
    // temporal separation. Its purpose it to limit the Q durting re-entry.  We
    // assume that this burn occurs at full-thrust, with AoA=0, so the only var-
    // iable param is the burn duration:
    // XXX: Do we always need FullThrust here?
    Pressure              m_entryBurnQ;
    Time                  m_entryBurnDur;

    // XXX: For the Final Descnet and Landing, we currently do NOT perform any
    // special maneuvers  to avoid the "ballistic target" point and fly to the
    // actual landing site. We just perform a continuous burn in order to land
    // with near-zero velocity. So we currently assume AoA=0 and a const  (but
    // subject to optimisation) "mu", given by the corresp ThrottlingLevel.  A
    // reasonable trigger is the FinalBurn altitude; the burn lasts until land-
    // ing (h=0), until v=0 (which means unsuccessful langing if "h" is above
    // the threshold), or until the propellant is exhausted:
    LenK                  m_finalBurnH;
    double                m_finalBurnThrtL;

    // So altogether: 8 to 11 params, depending on the number of BBBurnSinTheta
    // coeffs (from 1 to 4):
    constexpr static int  NP = 11;
    // [PropMassS,  CoastTime,    BBBurnDur,  BBBurnSinTheta[1..4],
    //  EntryBurnQ, EntryBurnDur, FinalBurnH, FinalBurnThrtL]

    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    RTLS1
    (
      // Stage Params:
      Mass           a_full_mass1,
      double         a_K1,                  // PropMass1 / FullMass1
      double         a_prop_rem1,
      Time           a_Isp_sl1,
      Time           a_Isp_vac1,
      ForceK         a_thrust_vac1,
      double         a_min_thrtl1,
      Len            a_diam,

      // Mission (Return-to-Launch-Site) Params:
      LenK           a_hS,
      LenK           a_lS,
      VelK           a_VrS,
      VelK           a_VhorS,

      // Integration and Output Params:
      Time           a_ode_integr_step,
      std::ostream*  a_os,
      int            a_log_level
    );

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
    Angle AoA(Time a_t, Angle a_psi) const;

    // "Thrust": Depends on the Mode, BurnRate and the Counter-Pressure:
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const;

    // "LVMass": Current Mass:
    Mass LVMass(StateV const& a_s, Time a_t) const;

    //-----------------------------------------------------------------------//
    // "NOMADEvaluator": Helper Class used in NOMAD Optimisation:            //
    //-----------------------------------------------------------------------//
    class        NOMADEvaluator;
    friend class NOMADEvaluator;

    //-----------------------------------------------------------------------//
    // "OPtRes" Struct: The result of "FindOptimalReturnCtls" below:         //
    //-----------------------------------------------------------------------//
    // It just contains the relevant flds taken from the main "RTLS1" class:
    //
    struct OptRes
    {
      // Data Flds:
      Mass       m_propMassS;
      Time       m_coastTime;
      Time       m_bbBurnDur;
      double     m_bbBurnSinTheta[NS];
      Pressure   m_entryBurnQ;
      Time       m_entryBurnDur;
      LenK       m_finalBurnH;
      double     m_finalBurnThrtL;

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
    FindOptimalAscentCtls
    (
      // All are are given via the ConfigFile.ini, as there are quite a few of
      // them:
      std::string const& a_config_ini,
      std::ostream*      a_os   // May be NULL
    );

  private:
    /*
    //=======================================================================//
    // "RunNOMAD":                                                           //
    //=======================================================================//
    // Helper invoked from "FindOptimalAscentCtls". Returns "true" on success
    // (then "a_init_vals" contains the optimal params),   "false" otherwise:
    //
    static bool RunNOMAD
    (
      // Main Optimisation Problem Setup:
      Ascent2 const*                      a_proto,
      std::array<bool,NP> const&          a_act_opts,
      std::vector<double>*                a_init_vals,
      std::vector<double> const&          a_lo_bounds,
      std::vector<double> const&          a_up_bounds,
      // Optimisation Constraints:
      int                                 a_max_evals,
      bool                                a_constr_q,
      bool                                a_constr_sep_q,
      bool                                a_constr_long_g,
      boost::property_tree::ptree const&  a_pt
    );
    */
  };
}
// End namespace SpaceBallistics
