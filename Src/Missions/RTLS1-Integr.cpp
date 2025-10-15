// vim:ts=2:et
//===========================================================================//
//                       "Missions/RTLS1-Integr.cpp":                        //
//               Return-to-Launch-Site Trajectory Integration                //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"

namespace SpaceBallistics
{
//===========================================================================//
// "RTLS1": Non-Default Ctor:                                                //
//===========================================================================//
// Initialises the Stage and Mission Params; the Ctl Params are set to the
// default vals:
RTLS1::RTLS1
(
  // Stage Params:
  Mass           a_max_full_mass1,
  Mass           a_empty_mass1,
  double         a_prop_rem1,
  Time           a_Isp_sl1,
  Time           a_Isp_vac1,
  ForceK         a_thrust_vac1,
  double         a_min_thrtl1,
  Len            a_diam,

  // Mission (Return-to-Launch-Site) Params:
  Mass           a_prop_massS,
  LenK           a_hS,
  LenK           a_lS,
  VelK           a_VrS,
  VelK           a_VhorS,

  // Integration and Output Params:
  Time           a_ode_integr_step,
  std::ostream*  a_os,
  int            a_log_level
)
: Base
  (
    // FullMassS (at Separation time):
    a_empty_mass1 + a_prop_massS,
    // The effective "K1":
    double(a_prop_massS  / (a_empty_mass1 + a_prop_massS)),
    // Re-calculate "PropRem":  Indeed, "a_prop_rem1" is based on the
    // FullPropMass, whereas we need one based on "a_prop_massS":
    double(a_prop_rem1 * (a_max_full_mass1 - a_empty_mass1) / a_prop_massS),
    a_Isp_sl1,
    a_Isp_vac1,
    a_thrust_vac1,
    a_min_thrtl1,
    a_diam,
    // Integration and Output Params:
    a_ode_integr_step,
    a_os,
    a_log_level
  ),
  // Mission Params (constant):
  m_hS            (a_hS),
  m_lS            (a_lS),
  m_VrS           (a_VrS),
  m_VhorS         (a_VhorS),
  // Ctl Params   (subject to Optimisation):
  m_propMassS     (a_prop_massS),
  // Other Ctl Params are set to their default vals as yet:
  m_coastTime     (),
  m_bbBurnDur     (),
  m_bbBurnSinTheta{ 0.0, 0.0, 0.0, 0.0 },
  m_entryBurnQ    (),
  m_entryBurnDur  (),
  m_finalBurnH    (),
  m_finalBurnThrtL(0.0)
{}

}
// End namespace SpaceBallistics
