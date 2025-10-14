// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Missions/RTLS1.h":                     //
//                Return-To-Launch-Site for Stage1 of an LV									 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Missions/LVBase.h"

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
  };
}
// End namespace SpaceBallistics
