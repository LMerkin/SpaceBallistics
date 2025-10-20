// vim:ts=2:et
//===========================================================================//
//                       "Missions/RTLS1-Integr.cpp":                        //
//           Return-To-(Launch/Landing)-Site Trajectory Integration          //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"
#include "LVBase.hpp"
#include "SpaceBallistics/PhysEffects/LVAeroDyn.hpp"
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
  Mass           a_max_full_mass1,
  double         a_full_k1,          // Based on "a_max_full_mass1"
  double         a_full_prop_rem1,   // ditto
  Time           a_Isp_sl1,
  Time           a_Isp_vac1,
  ForceK         a_thrust_vac1,
  double         a_min_thrtl1,
  Len            a_diam,

  // Mission (Return-To-(Launch/Landing)-Site) Params:
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
    // FullMassS (at Separation time):     EmptyMass1 + PropMassS:
    a_max_full_mass1 * (1.0 - a_full_k1) + a_prop_massS,

    // The new effective "K1":
    double(a_prop_massS /
          (a_max_full_mass1 * (1.0 - a_full_k1) + a_prop_massS)),

    // Re-calculate "PropRem":  Indeed, "a_prop_rem1" is based on the
    // FullPropMass, whereas we need one based on "a_prop_massS":
    double(a_full_prop_rem1 * a_max_full_mass1  * a_full_k1 / a_prop_massS),

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
  m_coastDur      (),
  m_bbBurnDur     (),
  m_bbBurnSinTheta{ 0.0, 0.0, 0.0, 0.0 },
  m_entryBurnQ    (),
  m_entryBurnDur  (),
  m_landBurnH     (),
  m_landBurnThrtL (0.0),
  // Transient Data:
  m_mode          (FlightMode::Coast),
  m_entryIgnTime  (NAN),
  m_landIgnTime   (NAN),
  m_finalTime     (NAN),
  m_eventStr      (),
  // NB: Constraints (subject to "max") must be initialisd to 0, not NAN:
  m_maxQ          (),
  m_sepQ          (),
  m_maxLongG      (0.0)
{
  // Checks:
  if (!(IsPos(m_hS) && IsPos(m_lS) && IsPos(m_VrS) && IsPos(m_VhorS)))
    throw std::invalid_argument("RTLS1::Ctor: Invalid Mission Param(s)");
}

//===========================================================================//
// "Run": Integrate the Return Trajectory:                                   //
//===========================================================================//
RTLS1::Base::RunRes RTLS1::Run()
{
  //-------------------------------------------------------------------------//
  // Run the integration FORWARDS from the Stage1 Separation point:          //
  //-------------------------------------------------------------------------//
  // (t=0 corresponds to Separation):
  LenK   rS     = R + m_hS;
  AngVel omegaS = 1.0_rad * m_VhorS / rS;
  assert(IsPos(omegaS));

  // The initial State Vector:
  Base::StateV s0 = std::make_tuple(rS, m_VrS, omegaS, 0.0_kg, 0.0_rad);

  // The RHS and the Call-Back Lambdas:
  auto rhs =
    [this](Base::StateV const&  a_s, Time a_t) -> Base::DStateV
    { return this->Base::ODERHS(a_s, a_t); };

  auto cb  =
    [this](Base::StateV*  a_s,  Time a_t, Time a_tau) -> bool
    { return this->ODECB (a_s,  a_t, a_tau); };

  // NB: Transient fields have been initialised in the Ctor. XXX: IMPORTANT:
  // It is currently assumed that "Run" is invoked only once per the object
  // lifetime!
  if (m_mode != FlightMode::Coast)
    throw std::runtime_error("RTLS1::Run: Repeated invocations not allowed");

  constexpr Time t0   = 0.0_sec;
  constexpr Time tMax = 3600.0_sec; // Certainly enough!
  Time           tEnd;  // Will be < tMax;
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
  }
  catch (Base::NearSingularityExn const& ns)
  {
    // We have reached a vicinity of the Singular Point:
    // This is a NORMAL (and moreover, a desirable) outcome. Compute a more
    // precise singular point position:
    return Base::LocateSingularPoint(ns);
  }
  // XXX: Any other exceptions are propagated to the top level, it's better
  // not to hide them...

  // Integration has run to completion, but not to the Singular Point:
  return Base::PostProcessRun(s0, tEnd);
}

//===========================================================================//
// ODE CallBack (invoked after the completion of each RKF5 step):            //
//===========================================================================//
// HERE FlightMode switching occurs, so this method is non-"const":
//
bool RTLS1::ODECB(StateV* a_s, Time a_t, Time a_tau)
{
  assert(a_s != nullptr && !IsPos(a_t));
  LenK   r             = std::get<0>(*a_s);
  LenK   h             = r - R;          // Altitude
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  Mass   spentPropMass = std::get<3>(*a_s);
  Angle  phi           = std::get<4>(*a_s);
  assert(IsPos(r) && !IsNeg(spentPropMass));

  //-------------------------------------------------------------------------//
  // Similar to "Base::ODERHS", check if we are approaching the singularity: //
  //-------------------------------------------------------------------------//
  VelK   Vhor          = r * omega / 1.0_rad;

  if (Vhor < Base::SingVhor)
  {
    omega              = AngVel(0.0);
    std::get<2>(*a_s)  = AngVel(0.0);
    Vhor               = VelK(0.0);
  }
  if (IsZero(omega) && Vr < Base::SingVr)
    // This is not an error -- just a stop integration now:
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

  // Constraints:
  double   M   = Base::Mach(atm, V);
  Density  rho = std::get<1>(atm);
  auto     V2  = To_Len_m(Sqr(V));
  Pressure Q   = 0.5 * rho * V2;

  m_maxQ       = std::max(m_maxQ, Q);
  assert(!IsNeg(longG));
  m_maxLongG   = std::max(m_maxLongG, longG);

  // Curr event for logging (if any):
  m_eventStr.clear();

  //-------------------------------------------------------------------------//
  // Mode Switching:                                                         //
  //-------------------------------------------------------------------------//
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

  // ExoAtmDesc -> EntryBurn: by Q:
  if (m_mode == FlightMode::ExoAtmDesc && Q >= m_entryBurnQ)
  {
    m_mode         = FlightMode::EntryBurn;
    m_entryIgnTime = a_t;
    m_eventStr     = "Exo-Atmospheric Descent Ends, Entry Burn Starts";
  }

  // EntryBurn -> EndoAtmDesc: by Time:
  if (m_mode == FlightMode::EntryBurn &&
      a_t    >= m_entryIgnTime + m_entryBurnDur)
  {
    m_mode     = FlightMode::EndoAtmDesc;
    m_eventStr = "Entry Burn Ends, Endo-Atmospheric Descent Starts";
  }

  // EndoAtmDesc -> LandBurn: by Altitude:
  if (m_mode == FlightMode::EndoAtmDesc && h <= m_landBurnH)
  {
    m_mode        = FlightMode::LandBurn;
    m_landIgnTime = a_t;
    m_eventStr    = "Endo-Atmospheric Descent Ends, Landing Burn Starts";
  }

  // Any State (but mainly Landing Burn) -> UNDEFINED: by SpentPropMass:
  if (spentPropMass >= Base::m_spendable1)
  {
    m_mode      = FlightMode::UNDEFINED;
    m_finalTime = a_t;
    m_eventStr  = "Stage1 Flamed-Out";
  }
  //-------------------------------------------------------------------------//
  // Integration Stopping Conds:                                             //
  //-------------------------------------------------------------------------//
  bool cont = true;

  // Do NOT continue if:
  // (*) we are on the surface (stopping by the Altitude);
  // (*) we are in the UNDEFINED mode (stopping by Ignition Time);
  // NB: both are NORMAL termination conds (as opposed to throwing
  //     the "StopNowExn"):
  if (!IsPos(h))
  {
    cont        = false;
    m_eventStr += "Integration Stopped @ H=0";
  }
  if (m_mode == FlightMode::UNDEFINED)
  {
    cont        = false;
    m_eventStr += "Integration Stopped @ MaxStartMass";
  }

  //-------------------------------------------------------------------------//
  // Log the Curr Event (if any):                                            //
  //-------------------------------------------------------------------------//
  // Angles in Degrees:
  Angle_deg psi_deg = To_Angle_deg(psi);
  Angle_deg aoa_deg = To_Angle_deg(aoa);

  if (!m_eventStr.empty() && Base::m_os != nullptr && Base::m_logLevel >= 2)
  {
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
  // Occurs with a 100 msec step, or if we are going to stop now:
  if (Base::m_os  != nullptr && Base::m_logLevel >= 3 &&
     (!cont || int(Round(double(a_t / 0.001_sec))) % 100 == 0))
  {
    // Thrust is more conveniently reported in kgf:
    auto tkg = thrust / g0K;

    *Base::m_os
      << a_t.Magnitude()     << '\t' << h.Magnitude()        << '\t'
      << L.Magnitude()       << '\t' << Vr.Magnitude()       << '\t'
      << Vhor.Magnitude()    << '\t' << V.Magnitude()        << '\t'
      << psi_deg.Magnitude() << '\t' << aoa_deg.Magnitude()  << '\t'
      << m.Magnitude()       << '\t' << ToString(m_mode)     << '\t'
      << tkg.Magnitude()     << '\t' << burnRate.Magnitude() << '\t'
      << Q.Magnitude()       << '\t' << M                    << '\t'
      << longG               << '\t' << a_tau.Magnitude()    << std::endl;
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
    // Both angles do not matter in this mode, so return 0s:
    return std::make_pair(0.0_rad, 0.0_rad);

  case FlightMode::BBBurn:
  {
    // We actually control THETA and derive the AoA from it:
    // sin(theta) is a cubic polynomial of "nt":
    double nt = double((a_t - m_coastDur) / m_bbBurnDur);

    // By construction, 0 <= nt <= 1, but the latter inequality may be broken
    // due to rounding errors etc, so:
    assert(nt >= 0.0);
    nt = std::min(nt, 1.0);

    double sinTheta =
      ((m_bbBurnSinTheta[3]  * nt + m_bbBurnSinTheta[2]) * nt +
        m_bbBurnSinTheta[1]) * nt + m_bbBurnSinTheta[0];
    sinTheta = std::min(1.0, std::max(-1.0, sinTheta));

    // We assume that during this burn, the thrust vector is always pointing
    // towards the "return", so cosTheta <= 0, thus:
    Angle theta(Pi<double> - ASin(sinTheta));  // theta in [Pi/2 .. 3*Pi/2]

    // FIXME: the AoA is not much relevant in this mode, since we disreagard
    // the aerodynamic forces during this Burn.
    // Define it formalluy as psi + AoA = theta:
    Angle aoa = theta - a_psi;
    return std::make_pair(aoa, theta);
  }

  case FlightMode::ExoAtmDesc:
  case FlightMode::EntryBurn:
  case FlightMode::EndoAtmDesc:
  case FlightMode::LandBurn:
  case FlightMode::UNDEFINED:
    // Assume we are descending "tail-first", so AoA = Pi and theta = Pi + psi;
    // XXX: there is currently no "AoA maneuver"  at the end  before the final
    // vertical landing:
    return std::make_pair(PI, PI + a_psi);

  default:
    assert(false);
    return std::make_pair(0.0_rad, 0.0_rad);
  }
  __builtin_unreachable();
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
RTLS1::AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const
{
  auto atm  = EAM::GetAtmConds(std::max(a_r - R, 0.0_km));

  // FIXME: In the "Coast" and "BBBurn" modes,  we completely disregard the
  // aerodynamic forces as yet -- in particular because there may be arbit-
  // rary AoAs for which we do not have a proper model yet; arguably, those
  // models are exo-atmospheric, so this assumption is not grossly imprecise:
  //
  if (m_mode == FlightMode::Coast || m_mode == FlightMode::BBBurn)
    return std::make_tuple(atm, ForceK(0.0), ForceK(0.0));

  // Otherwise: Endo-Atmospheric Movement:
  double M  = Base::Mach(atm, a_v);
  if (!IsFinite(M))
    return std::make_tuple(atm, ForceK(0.0), ForceK(0.0));

  // Generic Case:
  assert(M >= 0.0);

  // The Aerodynamic Force Main Term:
  ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * Base::m_crosS;

  // The Drag and Lift Coeffs, using the default model:
  // FIXME: There is no precise model for "tail-first" movement yet; we assume
  // that "cL" is the same as for the "head-first" movement,  whereas the "cD"
  // is just multiplied by 10 (???):
  double cD = LVAeroDyn::cD(M, a_AoA);
  double cL = LVAeroDyn::cL(M, a_AoA) * 10.0;

  return std::make_tuple(atm, cD * F, cL * F);
}

//===========================================================================//
// Propellant Burn Rate:                                                     //
//===========================================================================//
MassRate RTLS1::PropBurnRate(Time UNUSED_PARAM(a_t)) const
{
  switch (m_mode)
  {
  case FlightMode::Coast:
  case FlightMode::ExoAtmDesc:
  case FlightMode::EndoAtmDesc:
  case FlightMode::UNDEFINED:
    // Passive modes:
    return MassRate(0.0);

  case FlightMode::BBBurn:
  case FlightMode::EntryBurn:
    // XXX: Assume full MassRate here -- no throttling ctl yet:
    return Base::m_burnRateI1;

  case FlightMode::LandBurn:
    // XXX: Assume that only 1 engine is burning (of 9), at some throttled but
    // (as yet) constant rate:
    assert(0.0 <= m_landBurnThrtL  && m_landBurnThrtL <= 1.0);
    return Base::m_burnRateI1 / 9.0 * m_landBurnThrtL;

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
    assert(0.0 <= x && x <= 1.0);
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
  Mass spentPropMass = std::get<3>(a_s);
  assert(!IsNeg(spentPropMass));

  Mass m = Base::m_fullMass1 - spentPropMass;
  assert(IsPos(m));

  // In general, we must have m >= EmptyMass + UnSpendableMass, though this may
  // be broken in some boundary cases. Enforce this condition explicitly:
  m = std::max(m, Base::m_emptyMass1 + Base::m_unSpendable1);
  return m;
}
}
// End namespace SpaceBallistics
