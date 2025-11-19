// vim:ts=2:et
//===========================================================================//
//                       "Missions/RTLS1-Integr.cpp":                        //
//           Return-To-(Launch/Landing)-Site Trajectory Integration          //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include "LVBase.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"

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
  Mass           a_fpl_mass1,       // Assuming Full-Prop-Load
  double         a_fpl_k1,          // Based on "a_fpl_mass1"
  double         a_fpl_prop_rem1,   // ditto
  Time           a_Isp_sl1,
  Time           a_Isp_vac1,
  ForceK         a_thrust_vac1,
  double         a_min_thrtl1,
  Len            a_diam,

  // Mission (Return-To-(Launch/Landing)-Site) Params:
  Mass           a_prop_massS,
  LenK           a_hS,
  LenK           a_lS,
  VelK           a_VS,
  Angle          a_psiS,

  // Estimates and Constraints:
  VelK           a_dVhor_est,
  LenK           a_land_dl_limit,
  VelK           a_land_vel_limit,
  AccK           a_land_acc_limit,
  Pressure       a_Q_limit,
  double         a_longG_limit,
  bool           a_approx_land_burn,

  // Optimisation Ranges:
  double const   a_prop_massSN   [2],
  Time           a_min_coast_dur,
  double         a_bbb_theta_minN,
  double const   a_bbb_durN      [2],
  Time   const   a_entry_burn_dur[2],

  // Integration and Output Params:
  Time           a_ode_integr_step,
  std::ostream*  a_os,
  int            a_log_level
)
: Base
  (
    // FullMassS (at Separation time):     EmptyMass1 + PropMassS:
    a_fpl_mass1 * (1.0 - a_fpl_k1) + a_prop_massS,

    // The new effective "K1":
    double(a_prop_massS / (a_fpl_mass1 * (1.0 - a_fpl_k1) + a_prop_massS)),

    // Re-calculate "PropRem":  Indeed, "a_prop_rem1" is based on the
    // FullPropMass, whereas we need one based on "a_prop_massS":
    double(a_fpl_prop_rem1 * a_fpl_mass1 * a_fpl_k1 / a_prop_massS),

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
  m_hS                 (a_hS),
  m_lS                 (a_lS),
  m_VS                 (a_VS),
  m_psiS               (a_psiS),

  // Estimates and Limits:
  m_dVhorEst           (a_dVhor_est),
  m_landDLLimit        (a_land_dl_limit),
  m_landVelLimit       (a_land_vel_limit),
  m_landAccLimit       (a_land_acc_limit),
  m_QLimit             (a_Q_limit),
  m_longGLimit         (a_longG_limit),
  m_approxLandBurn     (a_approx_land_burn),

  // Memoised Params for a "Fully-Prop-Loaded" Stage1:
  m_fplMass1           (a_fpl_mass1),
  m_fplK1              (a_fpl_k1),
  m_fplPropRem1        (a_fpl_prop_rem1),
  m_diam               (a_diam),

  // Optimisation Ranges:
  m_propMassSRange     {a_prop_massSN   [0], a_prop_massSN   [1]},
  m_minCoastDur        (a_min_coast_dur),
  m_bbBurnThetaMinPi   (a_bbb_theta_minN),
  m_bbBurnDurRange     {a_bbb_durN      [0], a_bbb_durN      [1]},
  m_entryBurnDurRange  {a_entry_burn_dur[0], a_entry_burn_dur[1]},

  // Ctl Params are set to "NAN"s as yet:
  m_coastDur           (NAN),
  m_bbBurnDur          (NAN),
  m_bbBurnThrtL0       (NAN),
  m_bbBurnThrtL1       (NAN),
  m_bbBurnTheta0       (NAN),
  m_bbBurnTheta1       (NAN),
  m_entryBurnQ         (NAN),
  m_entryBurnDur       (NAN),
  m_entryBurnThrtL0    (NAN),
  m_entryBurnThrtL1    (NAN),
  m_landBurnH          (NAN),
  m_landBurnGamma      (NAN),
  // Transient Data:
  m_mode               (FlightMode::UNDEFINED),
  m_entryIgnTime       (NAN),
  m_eventStr           (),
  m_nextOutputTime     (),
  // NB: Constraints (subject to "max") must be initialised to 0, not NAN:
  m_maxQ               (),
  m_sepQ               (),
  m_maxLongG           (0.0)
{
  // Checks:
  if (!(IsPos(m_hS)     && IsPos(m_lS)       && IsPos(m_VS) && IsPos(m_psiS) &&
        m_psiS < PI_2   && IsPos(m_dVhorEst) && IsPos(m_landDLLimit)         &&
        IsPos(m_QLimit) && m_longGLimit > 0.0))
    throw std::invalid_argument("RTLS1::Ctor: Invalid Param(s)");

  // "m_propMass1" (in the Base) should be equal to "a_prop_massS" up to
  // rounding errors:
  assert (Base::m_propMass1.ApproxEquals(a_prop_massS));

  // Obviously, "a_propMassS" must be within the following limits:
  if (!(IsPos(a_prop_massS) && a_prop_massS < a_fpl_mass1 * a_fpl_k1))
    throw std::invalid_argument("RTLS1::Ctor: Invalid PropMassS");

  // XXX: The Optimisation Ranges are chcked at the point of use...
}

//===========================================================================//
// "Run": Integrate the Return Trajectory:                                   //
//===========================================================================//
RTLS1::Base::RunRes RTLS1::Run(FlightMode a_init_mode)
{
  // NB: The initial mode is either "Coast" (Main Run) or "EndoAtmDesc" (Calib-
  // ration Run):
  assert(a_init_mode == FlightMode::Coast ||
         a_init_mode == FlightMode::EndoAtmDesc);
  m_mode = a_init_mode;

  //-------------------------------------------------------------------------//
  // Run the integration FORWARDS from the Stage1 Separation point:          //
  //-------------------------------------------------------------------------//
  // (t=0 corresponds to Separation):
  LenK   rS     = R + m_hS;
  VelK   VrS    = m_VS * Sin(m_psiS);
  AngVel omegaS = m_VS * Cos(m_psiS) / rS * 1.0_rad;
  Angle  phiS   = m_lS / R * 1.0_rad;
  assert(IsPos(VrS) && IsPos(omegaS) && IsPos(phiS));

  // The initial State Vector:
  Base::StateV s0 = std::make_tuple(rS, VrS, omegaS, 0.0_kg, phiS);

  // The RHS and the Call-Back Lambdas (IsAscent=false):
  auto rhs =
    [this](Base::StateV const&  a_s, Time a_t) -> Base::DStateV
    { return this->Base::ODERHS(a_s, a_t, false); };

  auto cb  =
    [this]
    (
      Base::StateV*        a_s,
      Time                 a_t,
      Base::DStateV const& a_ds,
      Base::StateV  const& a_prev_s,
      Time                 a_prev_t,
      Base::DStateV const& a_prev_ds
    )
    -> bool
    { return this->ODECB (a_s, a_t, a_ds, a_prev_s, a_prev_t, a_prev_ds); };

  constexpr Time t0   = 0.0_sec;
  constexpr Time tMax = 3600.0_sec; // Certainly enough!
  Time           tEnd;              // Will be < tMax
  m_nextOutputTime    = 0.0_sec;    // Just to make sure

  Base::RunRes   res;               // "Error" as yet
  try
  {
    //-----------------------------------------------------------------------//
    // Actually Run the Integrator!                                          //
    //-----------------------------------------------------------------------//
    tEnd =
      RKF5(&s0,    t0, tMax, rhs,
           Base::m_odeIntegrStep,  Base::m_odeIntegrStep,
           Base::ODERelPrec,       &cb,  Base::m_os);
    assert(tEnd <= tMax);

    // If we got here: Integration has run to completion, but not to the Singu-
    // lar Point. Invoke the "PostProcessRun" (IsAscent=false):
    res  = Base::PostProcessRun(s0, tEnd, false);
  }
  catch (Base::NearSingularityExn const& ns)
  {
    // We have reached a vicinity of the Singular Point:
    // This is a NORMAL (and moreover, a desirable) outcome. Compute a more
    // precise singular point position;  IsAscent=false:
    res  = Base::LocateSingularPoint(ns, false);
  }
  catch (LandBurnApproxExn const& lba)
  {
    res = LandBurnApprox(lba);
  }
  // XXX: Any other exceptions are propagated to the top level, it's better not
  // to hide them...

  // The result:
  if (Base::m_os != nullptr && Base::m_logLevel >= 3)
#   pragma omp critical(Output)
    *Base::m_os
      << "# t=" << res.m_T.Magnitude() << " sec: Integration Done, RC="
      << Base::ToString(res.m_rc)      << std::endl;
  return res;
}

//===========================================================================//
// ODE CallBack (invoked after the completion of each RKF5 step):            //
//===========================================================================//
// HERE FlightMode switching occurs, so this method is non-"const";
// some args are unused yet:
//
bool RTLS1::ODECB
(
  Base::StateV*       a_s,      Time a_t,      Base::DStateV const& a_ds,
  Base::StateV const& a_prev_s, Time a_prev_t, Base::DStateV const& a_prev_ds
)
{
  assert(a_s != nullptr && !IsNeg(a_t));
  LenK   r             = std::get<0>(*a_s);
  LenK   h             = r - R;          // Altitude
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  VelK   Vhor          = r * omega / 1.0_rad;
  Mass   spentPropMass = std::max(std::get<3>(*a_s), 0.0_kg);
  Angle  phi           = std::get<4>(*a_s);
  assert(IsPos(r) && !IsNeg(spentPropMass));

  //-------------------------------------------------------------------------//
  // Similar to "Base::ODERHS", check if we are approaching the Singularity: //
  //-------------------------------------------------------------------------//
  if (Abs(Vhor) < Base::SingVhor)
  {
    omega              = AngVel(0.0);
    std::get<2>(*a_s)  = AngVel(0.0);
    Vhor               = VelK(0.0);
  }
  if (IsZero(omega) && Abs(Vr) < Base::SingV)
    // This is not an error -- just a stop integration now. XXX: In theory, V=0
    // does not always imply a singularity; but in practice, it does:
    //
    throw Base::NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  //-------------------------------------------------------------------------//
  // Generic Case:                                                           //
  //-------------------------------------------------------------------------//
  // Compute the Non-Gravitational Acceleration Components:
  Mass          m = LVMass(*a_s, a_t);  // OK to use the "old" Mode yet
  VelK          V;
  Angle         psi, aoa;
  EAM::AtmConds atm;
  MassRate      burnRate;
  ForceK        thrust;
  double        lvAxis[2];
  AccK          ngAcc [2];
  double        longG;

  NonGravForces
    (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
     lvAxis, ngAcc, &longG);

  LenK     L   = R * phi / 1.0_rad; // Down-Range Earth Distance, <= 0

  //-------------------------------------------------------------------------//
  // Constraints:                                                            //
  //-------------------------------------------------------------------------//
  double   M   = Base::Mach(atm, V);
  Pressure pa  = std::get<0>(atm);
  Density  rho = std::get<1>(atm);
  auto     V2  = To_Len_m(Sqr(V));
  Pressure Q   = 0.5 * rho * V2;

  m_maxQ       = std::max(m_maxQ, Q);
  assert(!IsNeg(longG));
  m_maxLongG   = std::max(m_maxLongG, longG);

  //-------------------------------------------------------------------------//
  // Mode Switching:                                                         //
  //-------------------------------------------------------------------------//
  // Curr event for logging (if any):
  m_eventStr.clear();
  bool cont = true;

  // At the beginning: Stage1 separation:
  if (IsZero(a_t))
  {
    assert(m_mode == FlightMode::Coast);
    m_sepQ     = Q;
    m_eventStr = "Stage1 Separation, Coast Starts";
  }

  // Coast -> BBBurn: by Time:
  if (m_mode == FlightMode::Coast  && a_t >= m_coastDur)
  {
    m_mode     = FlightMode::BBBurn;
    m_eventStr = "Coast Ends, Boost-Back Burn Starts";
  }
  
  // BBBurn -> ExoAtmDesc: by Time:
  if (m_mode == FlightMode::BBBurn && a_t >= m_coastDur + m_bbBurnDur)
  {
    m_mode     = FlightMode::ExoAtmDesc;
    m_eventStr = "Boost-Back Burn Ends, Exo-Atmospheric Descent Starts";
  }

  // ExoAtmDesc -> EntryBurn: by Q
  // (XXX: may not happen at all if "EntryBurnQ" is set incorrectly):
  if (m_mode == FlightMode::ExoAtmDesc && Q >= m_entryBurnQ)
  {
    m_mode         = FlightMode::EntryBurn;
    m_entryIgnTime = a_t;
    m_eventStr     = "Exo-Atmospheric Descent Ends, Entry Burn Starts";
  }

  // EntryBurn -> EndoAtmDesc: by Time
  // (if there was an "EntryBurn" in the 1st place):
  if (m_mode == FlightMode::EntryBurn &&
      a_t    >= m_entryIgnTime + m_entryBurnDur)
  {
    m_mode     = FlightMode::EndoAtmDesc;
    m_eventStr = "Entry Burn Ends, Endo-Atmospheric Descent Starts";
  }

  // In any mode (but generically in "EntryBurn"): Set the "LandBurn" params
  // which will (later or immediately) trigger the "LandBurn".
  // NB: These params would already be set in the Calibration mode, in which
  // case "SetLandBurnParams" is not invoked:
  //
  if (h <= MaxLandBurnH && !IsFinite(m_landBurnH))
  {
    // Then "m_landBurnGamma" is not set either:
    assert(!IsFinite(m_landBurnGamma));

    if (m_approxLandBurn)
      // If only APPROXIMATE LandBurn evaluation is required, immediately return
      // an estimate of the outcome, via an exception:
      throw LandBurnApproxExn(Vr, Vhor, spentPropMass, phi, a_t);

    // Otherwise: Set both of them using pre-calibrated functions:
    SetLandBurnParams(*a_s, a_t, a_ds, a_prev_s, a_prev_t, a_prev_ds);

    assert(!IsNeg(m_landBurnH)    && m_landBurnH     <= MaxLandBurnH &&
           0.0 <= m_landBurnGamma && m_landBurnGamma <= 1.0);
  }

  // EndoAtmDesc -> LandBurn: by Altitude;
  // but more generally, if we have missed the "EntryBurn" (because the MaxQ
  // was naturally below the threshold), then we would miss "EndoAtmDesc" as
  // well. So for robustness, we allow switching to "LandBurn" from ANY mode,
  // simply by the altitude:
  if (h <= m_landBurnH && m_mode != FlightMode::LandBurn)
  {
    m_mode     = FlightMode::LandBurn;
    m_eventStr = "Endo-Atmospheric Descent Ends, Landing Burn Starts";
  }

  // LandBurn -> UNDEFINED:
  // Switching Off the Engine when the altitude is less than the target one
  // (XXX: in addition,  the Engine(s) are switched off in any mode when we
  // run out of propellant -- eg in "BurnRate").
  // But more generally, we can do so in ANY mode (for robustness).  Normally,
  // we expect that integration should terminate via the "NearSingularityExn":
  //
  if (h <= 0.001_km && m_mode != FlightMode::UNDEFINED)
  {
    m_mode     = FlightMode::UNDEFINED;
    m_eventStr = "Stage1 Landed";
    cont       = false;       // Stop now!
  }

  //-------------------------------------------------------------------------//
  // Log the Curr Event (if any):                                            //
  //-------------------------------------------------------------------------//
  // Angles in Degrees:
  Angle_deg psi_deg = To_Angle_deg(psi);
  Angle_deg aoa_deg = To_Angle_deg(aoa);

  if (!m_eventStr.empty() && Base::m_os != nullptr && Base::m_logLevel >= 3)
  {
#   pragma omp critical(Output)
    *Base::m_os
      << "# t="           << a_t.Magnitude()
      << " sec: "         << m_eventStr
      << ": Mode="        << ToString(m_mode)
      << ", h="           << h.Magnitude()
      << " km, L="        << L.Magnitude()
      << " km, V="        << V.Magnitude()
      << " km/sec, Vr="   << Vr.Magnitude()
      << " km/sec, Vhor=" << Vhor.Magnitude()
      << " km/sec, psi="  << psi_deg.Magnitude()
      << " deg, AoA="     << aoa_deg.Magnitude()
      << " deg, LVMass="  << m.Magnitude()
      << " kg, Q="        << Q.Magnitude()
      << " Pa, M="        << M
      << ", LongG="       << longG
      << std::endl;
  }
  //-------------------------------------------------------------------------//
  // Main Output:                                                            //
  //-------------------------------------------------------------------------//
  // Occurs with a 0.1 sec step, or if we are going to stop now:
  if (Base::m_os != nullptr && Base::m_logLevel >= 4 &&
     (!cont || a_t >= m_nextOutputTime))
  {
    m_nextOutputTime += 0.1_sec;

    // Thrust is more conveniently reported in kgf:
    auto tkgf = thrust / g0K;

    // Time Step (for the step done):
    Time dt   = a_t - a_prev_t;
    assert(IsPos(dt));

#   pragma omp critical(Output)
    *Base::m_os
      << a_t.Magnitude()     << '\t' << h.Magnitude()        << '\t'
      << L.Magnitude()       << '\t' << Vr.Magnitude()       << '\t'
      << Vhor.Magnitude()    << '\t' << V.Magnitude()        << '\t'
      << psi_deg.Magnitude() << '\t' << aoa_deg.Magnitude()  << '\t'
      << m.Magnitude()       << '\t' << ToString(m_mode)     << '\t'
      << tkgf.Magnitude()    << '\t' << burnRate.Magnitude() << '\t'
      << Q.Magnitude()       << '\t' << M                    << '\t'
      << longG               << '\t' << pa.Magnitude()       << '\t'
      << dt.Magnitude()      << std::endl;
  }
  return cont;
}

//===========================================================================//
// Angle-of-Attack and Thrust Vector Elevation:                              //
//===========================================================================//
//        AoA   theta
std::pair<Angle,Angle> RTLS1::AoA(Time a_t, Angle a_psi) const
{
  switch (m_mode)
  {
  case FlightMode::Coast:
    // Both angles do not matter in this mode, so return AoA=0, theta=psi
    // (we still preserve the invariant theta = AoA + psi, that is,  here
    // we are still flying Head-First):
    return std::make_pair(0.0_rad, a_psi);

  case FlightMode::BBBurn:
  {
    // We actually control THETA and derive the AoA from it. It is a linear
    // function of the relative BBBurn time so far:
    double nt = double((a_t - m_coastDur) / m_bbBurnDur);
    nt = std::min(1.0, std::max(0.0, nt));

    Angle  theta = m_bbBurnTheta0 * (1.0 - nt) + m_bbBurnTheta1 * nt;
    // XXX: The following cond holds, but up to rounding errors:
    // MinBBBurnThetaPiFrac * PI <= theta <= PI
    theta = std::min(theta, PI);

    // XXX: the AoA is not much relevant in this mode, since we disreagard
    // the aerodynamic forces during this Burn; define it formally via the
    // invariant: psi + AoA = theta, w/o normalisation:
    Angle aoa = theta - a_psi;
    return std::make_pair(aoa, theta);
  }

  case FlightMode::ExoAtmDesc:
  case FlightMode::EntryBurn:
  case FlightMode::EndoAtmDesc:
  case FlightMode::LandBurn:
  case FlightMode::UNDEFINED:
    // Assume we are descending Tail-First, so AoA = Pi  and theta = Pi + psi;
    // XXX: there is currently no "overfly maneuver" at the end before the fi-
    // nal vertical landing:
    return std::make_pair(PI, PI + a_psi);

  default:
    assert(false);
    return std::make_pair(0.0_rad, a_psi);
  }
  __builtin_unreachable();
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
RTLS1::AeroDynForces(LenK a_r, VelK a_v, Angle DEBUG_ONLY(a_AoA)) const
{
  auto atm  = EAM::GetAtmConds(std::max(a_r - R, 0.0_km));

  // NB: In the "Coast" and "BBBurn" modes, we completely disregard the aerody-
  // namic forces as yet -- in particular because there may be arbitrary AoAs
  // for which we do not have a proper model yet; arguably,   those modes are
  // exo-atmospheric, so this assumption is reasonably precise:
  //
  if (m_mode == FlightMode::Coast || m_mode == FlightMode::BBBurn)
    return std::make_tuple(atm, ForceK(0.0), ForceK(0.0));

  // Otherwise: Endo-Atmospheric Movement. XXX: Here we must have AoA=Pi for
  // now:
  assert(a_AoA == PI);

  double M  = Base::Mach(atm, a_v);
  assert(IsFinite(M) && M >= 0.0);

  // The following is a rough but reasonable approximation for the "cD" during
  // the "tail-first" motion:
  double cD =
    1.40 + 0.25 * Exp(-Sqr(M-1.0) / 0.12) - 0.15 * TanH(0.8 * (M - 1.4));
  assert(cD > 0.0);

  // The Aerodynamic Drag Force:
  ForceK FD =
    cD * 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * Base::m_crosS;

  // The Lift Force @ AoA=0 is assumed to be 0. Thus:
  return std::make_tuple(atm, FD, ForceK(0.0));
}

//===========================================================================//
// Propellant Burn Rate:                                                     //
//===========================================================================//
MassRate RTLS1::PropBurnRate(Mass a_m, LenK a_h, Time a_t) const
{
  // IMPORTANT: Irrespective to the Mode, if we are out of Propellant, the
  // BurnRate is obviously 0:
  if (a_m <= Base::m_emptyMass1 + Base::m_unSpendable1)
    return MassRate(0.0);

  // Otherwise:
  switch (m_mode)
  {
  case FlightMode::Coast:
  case FlightMode::ExoAtmDesc:
  case FlightMode::EndoAtmDesc:
  case FlightMode::UNDEFINED:
    // Passive modes:
    return MassRate(0.0);

  case FlightMode::BBBurn:
  {
    // Linear function of the relative BBBurn time so far:
    double   nt  = double((a_t - m_coastDur) / m_bbBurnDur);
    nt           = std::min(1.0, std::max(0.0, nt));
    assert(0.0 <= m_bbBurnThrtL1 && m_bbBurnThrtL1 <= m_bbBurnThrtL0);
    return
      Base::m_burnRateI1 * BBBurnEngPart *
      (m_bbBurnThrtL0 * (1.0 - nt) + m_bbBurnThrtL1 * nt);
  }

  case FlightMode::EntryBurn:
  {
    // Linear function of the relative EntryBurn time so far:
    double nt = double((a_t - m_entryIgnTime) / m_entryBurnDur);
    nt        = std::min(1.0, std::max(0.0, nt));
    assert(0.0 <= m_entryBurnThrtL1 && m_entryBurnThrtL1 <= m_entryBurnThrtL0);
    return
      Base::m_burnRateI1 * EntryBurnEngPart *
      (m_entryBurnThrtL0 * (1.0 - nt) + m_entryBurnThrtL1 * nt);
  }

  case FlightMode::LandBurn:
  {
    // XXX: In this case, the BurnRate is an up-convex function of the Altitude,
    // to make sure that we begin the LandBurn with 100% throttling level  (but
    // still with "LandBurnEngPart"!) for sufficiently quick deceleration,  and
    // then reduce the throttling level to "m_minThrtL1" for the smallest decel-
    // eration at touch-down:
    if (IsPos(m_landBurnH))
    {
      // Generic Case:
      double x     =
        std::pow(std::max(0.0, std::min(1.0, double(a_h / m_landBurnH))),
                 m_landBurnGamma);
      double thrtL = Base::m_minThrtL1 * (1.0 - x) + x;

      return Base::m_burnRateI1 * LandBurnEngPart * thrtL;
    }
    else
      // There is no LandBurn:
      return MassRate(0.0);
  }

  default:
    assert(false);
    return MassRate(0.0);
  }
  __builtin_unreachable();
}

//===========================================================================//
// Thrust:                                                                   //
//===========================================================================//
ForceK RTLS1::Thrust(MassRate a_burn_rate, Pressure a_p) const
{
  assert(!(IsNeg(a_burn_rate) || IsNeg(a_p)));

  // For optimisation: If BurnRate is 0 (in any Mode), then so is Thrust:
  if (IsZero(a_burn_rate))
    return ForceK(0.0);

  switch (m_mode)
  {
  case FlightMode::Coast:
  case FlightMode::ExoAtmDesc:
  case FlightMode::EndoAtmDesc:
  case FlightMode::UNDEFINED:
    // These modes are passive:
    assert(IsZero(a_burn_rate));
    return ForceK(0.0);

  case FlightMode::BBBurn:
    // XXX: This Burn is considered to be Exo-Atmospheric,   so we currently
    // assume there is no static counter-pressure, and the Isp is the Vacuum
    // one:
    return a_burn_rate * Base::m_IspVac1 * g0K;

  case FlightMode::EntryBurn:
  case FlightMode::LandBurn:
  {
    // Here the static counter-pressure is taken into account:
    double   x   = double(a_p / EAM::P0);
    assert(0.0  <= x &&  x <= 1.0);
    Time     Isp = x   * Base::m_IspSL1 + (1.0 - x) * Base::m_IspVac1;
    return a_burn_rate * Isp * g0K;
  }

  default:
    assert(false);
    return ForceK(0.0);
  }
  __builtin_unreachable();
}

//===========================================================================//
// Current LV (actually Stage1) Mass:                                        //
//===========================================================================//
Mass RTLS1::LVMass(Base::StateV const& a_s, Time UNUSED_PARAM(a_t)) const
{
  Mass spentPropMass = std::max(std::get<3>(a_s), 0.0_kg);
  assert(!IsNeg(spentPropMass));

  Mass m = Base::m_fullMass1 - spentPropMass;
  assert(IsPos(m));

  // In general, we must have m >= EmptyMass + UnSpendableMass, though this may
  // be broken in some edge cases due to rounding errors:
  m = std::max(m, Base::m_emptyMass1 + Base::m_unSpendable1);
  return m;
}
}
// End namespace SpaceBallistics
