// vim:ts=2:et
//===========================================================================//
//                       "Missions/RTLS1-Integr.cpp":                        //
//           Return-to-(Launch/Landing)-Site Trajectory Integration          //
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
  Mass           a_empty_mass1,
  double         a_prop_rem1,
  Time           a_Isp_sl1,
  Time           a_Isp_vac1,
  ForceK         a_thrust_vac1,
  double         a_min_thrtl1,
  Len            a_diam,

  // Mission (Return-to-(Launch/Landing)-Site) Params:
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
  m_coastDur      (),
  m_bbBurnDur     (),
  m_bbBurnSinTheta{ 0.0, 0.0, 0.0, 0.0 },
  m_entryBurnQ    (),
  m_entryBurnDur  (),
  m_finalBurnH    (),
  m_finalBurnThrtL(0.0),
  // Transient Data:
  m_mode          (FlightMode::Coast),
  m_entryIgnTime  (NAN),
  m_finalIgnTime  (NAN),
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
  LenK   h             = r - R;             // Altitude
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
  Mass          m      = LVMass(*a_s, a_t);  // Using the "old" Mode yet!
  VelK          V;
  Angle         psi, aoa;
  EAM::AtmConds atm;
  MassRate      burnRate;
  ForceK        thrust;
  double        lvAxis[2];
  AccK          ngAcc [2];
  NonGravForces
    (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
     lvAxis, ngAcc);

  LenK      L          = R * phi / 1.0_rad; // Down-Range Earth Distance, <= 0
  Angle_deg psi_deg    = To_Angle_deg(psi);
  Angle_deg aoa_deg    = To_Angle_deg(aoa);

  // AeroDynamic Conditions:
  double   M           = Base::Mach(atm, V);
  Density  rho         = std::get<1>(atm);
  auto     V2          = To_Len_m(Sqr(V));
  Pressure Q           = 0.5 * rho * V2;

  // Output at the beginning:
  if (IsZero(a_t) && Base::m_os != nullptr && m_logLevel >= 1)
  {
    assert(m_mode == FlightMode::Coast);
    *Base::m_os << "# t=0 sec, h=" << h.Magnitude()
                << " km, V="       << V.Magnitude() << " km/sec, Mass="
                << m.Magnitude()   << " kg"         << std::endl;
  }

  //-------------------------------------------------------------------------//
  // Switching Coast -> BBBurn:                                              //
  //-------------------------------------------------------------------------//
  if (m_mode == FlightMode::Coast && a_t >= m_coastDur)
  {
    m_mode = FlightMode::BBBurn;

    if (Base::m_os != nullptr && Base::m_logLevel >= 2)
      *Base::m_os
        << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
        << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
        << " km/sec, psi="              << psi_deg.Magnitude()
        << " deg, m=" << m.Magnitude()  << " kg: "
           "Boost-Back Burn Starts"     << std::endl;
  }

/*
      UNDEFINED = 0,
      Coast     = 1,  // From separation to BoostBackBurn (BBBurn): passive
      BBBurn    = 2,  //
      ExAtmDesc = 3,  // Exo-Atmospheric  Descent:                  passive
      EntryBurn = 4,  // Slowing-Down before Re-Entering the Atmosphere
      EnAtmDesc = 5,  // Endo-Atmospheric Descent:                  passive
      FinalBurn = 6   // 
*/
  return true;
}
}
// End namespace SpaceBallistics
