// vim:ts=2:et
//===========================================================================//
//                    "Src/Missions/Ascent2-Integr.cpp":                     //
//     Ascent-to-Orbit for a "Model" 2-Stage LV: Trajectory Integration      //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"
#include "LVBase.hpp"
#include "SpaceBallistics/PhysEffects/LVAeroDyn.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"

namespace SpaceBallistics
{
//===========================================================================//
// "Ascent2": Non-Default Ctor:                                              //
//===========================================================================//
// Initialises the LV and Mission Params; the Ctl Params are set to the default
// vals:
Ascent2::Ascent2
(
  // Stage2:
  double          a_K2,
  double          a_prop_rem2,
  Time            a_Isp_vac2,
  ForceK          a_thrust_vac2,
  double          a_min_thrtl2,
  Angle_deg       a_max_aoa2,

  // Stage1:
  double          a_K1,
  double          a_prop_rem1,
  Time            a_Isp_sl1,
  Time            a_Isp_vac1,
  ForceK          a_thrust_vac1,
  double          a_min_thrtl1,
  Angle_deg       a_max_aoa1,

  // Over-All:
  double          a_alpha1,
  Mass            a_max_start_mass,
  Mass            a_fairing_mass,
  Len             a_diam,
  Mass            a_payload_mass,

  // Mission Params:
  LenK            a_h_perigee,
  LenK            a_h_apogee,
  Angle_deg       a_incl,
  Angle_deg       a_launch_lat,

  // Integration/Output Params:
  Time            a_ode_integr_step,
  std::ostream*   a_os,
  int             a_log_level
)
: Base
  (
    // FullMass1:
    (a_max_start_mass  - a_fairing_mass - a_payload_mass) / (1.0 + a_alpha1) *
     a_alpha1,
    // Other Stage1 Params:
    a_K1,
    a_prop_rem1,
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
  //-------------------------------------------------------------------------//
  // Over-All and some Stage1 Params:                                        //
  //-------------------------------------------------------------------------//
  // Const Over-All Params (not affected by Optimisation):
  m_maxStartMass  (a_max_start_mass),
  m_fairingMass   (a_fairing_mass),

  // Non-Const Over-All Params (may be updated in the course of Optimisation):
  m_alpha1        (a_alpha1),
  m_payLoadMass   (a_payload_mass),

  //-------------------------------------------------------------------------//
  // Stage2:                                                                 //
  //-------------------------------------------------------------------------//
  // Const Stage2 Params (not affected by Optimisation):
  m_K2            (a_K2),
  m_propRem2      (a_prop_rem2),
  m_IspVac2       (a_Isp_vac2),
  m_minThrtL2     (a_min_thrtl2),

  // Non-Const Stage2 Params (may be updated in the course of Optimisation):
  m_fullMass2     ((m_maxStartMass  - m_fairingMass - m_payLoadMass) /
                   (1.0 + m_alpha1)),
  m_emptyMass2    (m_fullMass2   * (1.0 - m_K2)),
  m_propMass2     (m_fullMass2   * m_K2),
  m_unSpendable2  (m_propMass2   * m_propRem2),
  m_spendable2    (m_propMass2   - m_unSpendable2),
  m_thrustVacI2   (a_thrust_vac2),
  m_burnRateI2    (m_thrustVacI2 / (m_IspVac2 * g0K)),
  m_T2            (m_spendable2  / m_burnRateI2),

  //-------------------------------------------------------------------------//
  // Mission Params: Ininitiased in the Ctor Body:                           //
  //-------------------------------------------------------------------------//
  m_Rins          (0.0_km),
  m_Vins          (0.0),

  //-------------------------------------------------------------------------//
  // Flight Ctl Params:                                                      //
  //-------------------------------------------------------------------------//
  // Stage2:
  m_thrustMult2   (1.0),
  m_aMu2          (0.0),
  m_bMu2          (0.0),
  m_bHat2         (0.0),
  m_muHat2        (1.0),

  m_maxAoA2       (To_Angle(a_max_aoa2)),
  m_aAoA2         (0.0),
  m_bAoA2         (0.0),
  m_aAoAHat2      (0.0),
  m_bAoAHat2      (0.0),

  // Ballistic Gap:
  m_TGap          (0.0_sec),

  // Stage1:
  m_thrustMult1   (1.0),
  m_aMu1          (0.0),
  m_bMu1          (0.0),
  m_bHat1         (0.0),
  m_muHat1        (1.0),

  m_maxAoA1       (To_Angle(a_max_aoa1)),
  m_aAoA1         (0.0),
  m_bAoA1         (0.0),
  m_aAoAHat1      (0.0),
  m_bAoAHat1      (0.0),

  //-------------------------------------------------------------------------//
  // Transient Data:                                                         //
  //-------------------------------------------------------------------------//
  m_mode          (FlightMode::Burn2),
  m_ignTime2      (NAN),
  m_fairingSepTime(NAN),
  m_cutOffTime1   (NAN),
  m_ignTime1      (NAN),
  m_eventStr      (),
  m_nextOutputTime(),
  m_maxQ          (0.0),
  m_sepQ          (0.0),
  m_maxLongG      (0.0)
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  if (!(0.0 < m_K2            && m_K2 < 1.0            &&
        IsPos(m_fullMass2)    && IsPos(m_emptyMass2)   &&
        IsPos(m_propMass2)    && IsPos(m_unSpendable2) &&
        IsPos(m_spendable2)   && IsPos(m_IspVac2)      &&
        IsPos(m_burnRateI2)   && IsPos(m_thrustVacI2)  &&
        IsPos(m_T2)           &&
        //
        0.0 <= m_minThrtL2    && m_minThrtL2  <= 1.0   &&
        IsPos(m_maxStartMass) && IsPos(m_fairingMass)  &&
        m_thrustMult2 > 0.0   && m_thrustMult1 > 0.0   &&
        !IsNeg(m_maxAoA2)     && !IsNeg(m_maxAoA1) ))
    throw std::invalid_argument("Ascent2::Ctor: Invalid LV Param(s)");

  // XXX: We assume that launch is NOT from the North or South Pole:
  if (!(a_h_perigee       >=  100.0_km && a_h_perigee <= a_h_apogee &&
        a_incl            >= 0.0_deg   && a_incl      <= 180.0_deg  &&
        Abs(a_launch_lat)  < 90.0_deg  && !IsNeg(m_payLoadMass)))
    throw std::invalid_argument("Ascent2::Ctor: Invalid Mission Param(s)");

  //-------------------------------------------------------------------------//
  // Mission Params:                                                         //
  //-------------------------------------------------------------------------//
  // Semi-Major Axis and the Perigee radius-vector:
  LenK a = R + 0.5 * (a_h_perigee + a_h_apogee);
  m_Rins = R +        a_h_perigee;

  // The required orbital velocity at the orbital insertion point (which is
  // always assumed to be the Perigee):
  VelK Vq = SqRt(K * (2.0 / m_Rins - 1.0 / a));

  // The Orbit Inclination and the Launch Site Latitude:
  double cosI   = Cos(To_Angle_rad(a_incl));
  double cosLat = Cos(To_Angle_rad(a_launch_lat));
  assert(cosLat > 0.0);

  if (std::fabs(cosI) > cosLat)
    throw std::invalid_argument("Ascent2::Ctor: Inclination Unreachable");

  // "A" is the Launch Azimuth:
  double sinA   = cosI / cosLat;
  assert(std::fabs(sinA) <= 1.0);

  // We will assume -Pi/2 <= A <= Pi/2, so cosA >= 0:
  double cosA   = SqRt(1.0 - Sqr(sinA));

  // The Earth Rotation Velocity at the Launch Latitude:
  VelK   ERV    =
      TwoPi<double> * R / BodyData<Body::Earth>::SiderealRotationPeriod *
      cosLat;

  // The components of the velocity (in the Launch COS) which must be attained
  // by the LV at the orbital insertion point, taking the Earth rotation into
  // account:
  VelK   Vlat   = Vq * sinA - ERV;
  VelK   Vlong  = Vq * cosA;
  // So the LV velocity at the orbital insertion point must be:
  m_Vins        = SqRt(Sqr(Vlat) + Sqr(Vlong));
}

//===========================================================================//
// "Run": Integrate the Ascent Trajectory:                                   //
//===========================================================================//
Ascent2::Base::RunRes Ascent2::Run()
{
  //-------------------------------------------------------------------------//
  // Run the integration BACKWARDS from the Orbital Insertion point          //
  //-------------------------------------------------------------------------//
  // (@ t=0):
  // Angular LV velocity at the orbital insertion point:
  AngVel omega0 = 1.0_rad * m_Vins / m_Rins;

  // The initial Radial Velocity is 0, and the final delta of Spent Mass
  // is also 0, and the initial Polar Angle is 0:
  Base::StateV s0 = std::make_tuple(m_Rins, VelK(0.0), omega0, 0.0_kg, 0.0_rad);

  // The RHS and the Call-Back Lambdas (IsAscent=true):
  auto rhs =
    [this](Base::StateV const&  a_s, Time a_t) -> Base::DStateV
    { return this->Base::ODERHS(a_s, a_t, true); };

  auto cb  =
    [this](Base::StateV*  a_s,  Time a_t, Time a_dt) -> bool
    { return this->ODECB (a_s,  a_t, a_dt); };

  // NB: Transient fields have been initialised in the Ctor. XXX: IMPORTANT:
  // It is currently assumed that "Run" is invoked only once per the object
  // lifetime!
  if (m_mode != FlightMode::Burn2)
    throw std::runtime_error("Ascent2::Run Repeated invocations not allowed");

  // The ascent maximum duration (which is certainly enough) is 1 hour:
  constexpr Time t0   = 0.0_sec;
  constexpr Time tMin = -3600.0_sec;  // Certainly enough!
  Time           tEnd;                // Will be > tMin
  m_nextOutputTime    = 0.0_sec;      // Just to make sure

  Base::RunRes   res;                 // "Error" as yet
  try
  {
    //-----------------------------------------------------------------------//
    // Actually Run the Integrator!                                          //
    //-----------------------------------------------------------------------//
    tEnd =
      RKF5(&s0,    t0, tMin, rhs,
           -Base::m_odeIntegrStep, -Base::m_odeIntegrStep,
            Base::ODERelPrec, &cb,  Base::m_os);
    assert(tEnd >= tMin);

    // If we got here: Integration has run to completion, but not to the Singu-
    // lar Point. Invoke the "PostProcessRun" (IsAscent=false):
    res  = Base::PostProcessRun(s0, tEnd, false);
  }
  catch (Base::NearSingularityExn const& ns)
  {
    // We have reached a vicinity of the Singular Point:
    // This is a NORMAL (and moreover, a desirable) outcome. Compute a more
    // precise singular point position;  IsAscent=true:
    res  = Base::LocateSingularPoint(ns, true);
  }
  catch (FallingBackExn const& fb)
  {
    // No point in continuing: we are falling back to Earth instead of ascend-
    // ing to orbit:
    VelK V = SqRt(Sqr(fb.m_Vr) + Sqr(fb.m_Vhor));

    if (Base::m_os != nullptr)
#     pragma omp critical(Output)
      (*m_os)
        << "# ERROR: Falling Back to Earth: t="            << fb.m_t.Magnitude()
        << " sec,  h="         << (fb.m_r - R).Magnitude() << " km, l="
        << (R * fb.m_phi / 1.0_rad).Magnitude()            << " km, Vr="
        << fb.m_Vr.Magnitude() << " m/sec, V="             << V.Magnitude()
        << " km/sec, mode="    << ToString(m_mode)         << std::endl;

    // XXX: We do not propagate vals carried by "fb" to "RunRes" which remains
    // an error...
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
// HERE FlightMode switching occurs, so this method is non-"const":
//
bool Ascent2::ODECB(Base::StateV* a_s, Time a_t, Time a_dt)
{
  assert(a_s != nullptr && !IsPos(a_t));
  LenK   r             = std::get<0>(*a_s);
  LenK   h             = r - R;             // Altitude
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  VelK   Vhor          = r * omega / 1.0_rad;
  Mass   spentPropMass = std::get<3>(*a_s);
  Angle  phi           = std::get<4>(*a_s);

  // NB: "spentPropMass" decreases over time, so increases in Bwd time.
  // It is 0 at Orbit Insertion Time (t=0):
  assert(IsPos(r) && !IsNeg(spentPropMass));

  //-------------------------------------------------------------------------//
  // Check if we are approaching the Singularity:                            //
  //-------------------------------------------------------------------------//
  // XXX: Unlike Base::ODERHS, here we set the angular velocity to 0 if "Vhor"
  // (not |Vhor|!) is below the threshold, because we KNOW that this may only
  // happen near the Singular Point:
  //
  if (Vhor < Base::SingVhor)
  {
    omega              = AngVel();
    std::get<2>(*a_s)  = AngVel();
    Vhor               = VelK  ();
  }
  if (IsZero(omega) && Abs(Vr) < Base::SingV)
    // This is not an error -- just stop integration now. XXX: Small negative
    // vals of "Vr" are handled in this case as well. In theory, V=0 does not
    // always imply a singularity; but in practice, it does:
    //
    throw Base::NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  // Otherwise, we allow Vr < 0 ONLY during "Burn2" and if the Ascent trajectory
  // is above the orbital insertion point:  such Ascent profiles can indeed be
  // In all other cases, it means that we are falling back to Earth, so stop in-
  // tegration immediately:
  if (IsNeg(Vr) && !(m_mode == FlightMode::Burn2 && r > m_Rins))
  {
    // If it is small (< 10 m/sec), just correct it:
    if (Abs(Vr) < VelK(0.01))
      Vr = VelK(0.0);
    else
      // Otherwise, throw an exception; when handled it will invalidate this
      // solution:
      throw FallingBackExn(r, Vr, Vhor, spentPropMass, phi, a_t);
  }

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
  double        longG;

  Base::NonGravForces
    (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
     lvAxis, ngAcc, &longG);

  LenK     L   = R * phi / 1.0_rad; // Down-Range Earth Distance, <= 0
  // AeroDynamic Conditions:
  double   M   = Base::Mach(atm, V);
  Pressure pa  = std::get<0>(atm);
  Density  rho = std::get<1>(atm);
  auto     V2  = To_Len_m(Sqr(V));
  Pressure Q   = 0.5 * rho * V2;

  // Curr event for logging (if any):
  m_eventStr.clear();

  if (IsZero(a_t))
    m_eventStr = "Stage2 Cut-Off, Orbital Insertion";

  //-------------------------------------------------------------------------//
  // Switching Burn2 -> Gap:                                                 //
  //-------------------------------------------------------------------------//
  // Occurs according to "spentPropMass": when all Stage2 propellants (except
  // the Remnants) is spent:
  //
  if (m_mode == FlightMode::Burn2 && spentPropMass >= m_spendable2)
  {
    assert(!IsFinite(m_ignTime2));
    m_mode        = FlightMode::Gap;
    m_eventStr    = "Ballistic Gap Ends, Stage2 Ignition";
    m_ignTime2    = a_t;
    assert(IsNeg(m_ignTime2));

    // Furthermore, since the duration of the Ballistic Gap is known, at
    // this point we already know the Stage1 Cut-Off time:
    // NB:  HERE "m_TGap" is used:
    assert(!IsNeg(m_TGap) && !IsFinite(m_cutOffTime1));
    m_cutOffTime1 = a_t - m_TGap;
    assert(m_cutOffTime1 <= m_ignTime2);
  }

  //-------------------------------------------------------------------------//
  // Switching Gap -> Burn1:                                                 //
  //-------------------------------------------------------------------------//
  // Occurs merely by time delay (if the Gap is 0, this will happen immediately
  // after the above change Burn2 -> Gap):
  //
  if (m_mode  == FlightMode::Gap)
  {
    assert(IsNeg(m_cutOffTime1));
    if (a_t  <=  m_cutOffTime1)
    {
      m_mode     = FlightMode::Burn1;
      m_eventStr = "Stage1 Cut-Off and Separation, Ballistic Gap Starts";

      // XXX: We assume that Stage1 separation occures at the cut-off moment,
      // so record the corresp Q:
      m_sepQ   = Q;
    }
  }

  //-------------------------------------------------------------------------//
  // Switching Burn1 -> UNDEFINED:                                           //
  //-------------------------------------------------------------------------//
  // When all Stage1 propellants (except the Remnants) is spent:
  //
  if (m_mode == FlightMode::Burn1  &&
      spentPropMass >= m_spendable2 + Base::m_spendable1)
  {
    m_mode     = FlightMode::UNDEFINED;
    m_eventStr = "Stage1 Ignition, Lift-Off";
    assert(!IsFinite(m_ignTime1));
    m_ignTime1 = a_t;
    assert(m_ignTime1 < m_cutOffTime1 && m_cutOffTime1 <= m_ignTime2 &&
           IsNeg(m_ignTime2));
  }

  //-------------------------------------------------------------------------//
  // In any mode, detect the Fairing Separation Condition:                   //
  //-------------------------------------------------------------------------//
  // In the reverse time, it's when the Dynamic Pressure becomes HIGHER
  // that the threshold:
  if (!IsFinite(m_fairingSepTime) && Q >= FairingSepCond)
  {
    m_fairingSepTime = a_t;
    m_eventStr       = "Fairing Separation";
  }

  //-------------------------------------------------------------------------//
  // Monitor the Constraints:                                                //
  //-------------------------------------------------------------------------//
  // If there was any event, re-compute the LV Mass, Forces and Accelerations:
  if (!m_eventStr.empty())
  {
    m    = LVMass(*a_s, a_t);
    Base::NonGravForces
      (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
      lvAxis, ngAcc, &longG);
  }
  // Dynamic Pressure:
  m_maxQ     = std::max(Q, m_maxQ);

  // LongG:
  assert(!IsNeg(longG));
  m_maxLongG = std::max(longG,  m_maxLongG);

  //-------------------------------------------------------------------------//
  // Integration Stopping Conds:                                             //
  //-------------------------------------------------------------------------//
  bool cont = true;

  // Do NOT continue if:
  // (*) we are on the surface (stopping by "h");
  // (*) we are in the UNDEFINED mode (stopping by Ignition Time);
  // NB: both are NORMAL termination conds (as opposed to throwing
  //     the "StopNowExn"):
  if (!IsPos(h))
  {
    cont        = false;
    m_eventStr += ": Integration Stopped @ H <= 0";
  }
  if (m_mode == FlightMode::UNDEFINED)
  {
    cont        = false;
    m_eventStr += ": Integration Stopped @ MaxStartMass";
  }

  //-------------------------------------------------------------------------//
  // Log the Curr Event (if any):                                            //
  //-------------------------------------------------------------------------//
  // Angles in Degrees:
  Angle_deg psi_deg = To_Angle_deg(psi);
  Angle_deg aoa_deg = To_Angle_deg(aoa);

  if (!m_eventStr.empty() && Base::m_os != nullptr && Base::m_logLevel >= 3)
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

  //-------------------------------------------------------------------------//
  // Main Output:                                                            //
  //-------------------------------------------------------------------------//
  // Occurs with a 0.1 msec step, or if we are going to stop now:
  if (Base::m_os != nullptr && Base::m_logLevel >= 4 &&
     (!cont || a_t <= m_nextOutputTime))
  {
    m_nextOutputTime -= 0.1_sec;   // NB: Bwd Integration!

    // Thrust is more conveniently reported in kgf:
    auto tkg = thrust / g0K;

#   pragma omp critical(Output)
    *Base::m_os
      << a_t.Magnitude()     << '\t' << h.Magnitude()        << '\t'
      << L.Magnitude()       << '\t' << Vr.Magnitude()       << '\t'
      << Vhor.Magnitude()    << '\t' << V.Magnitude()        << '\t'
      << psi_deg.Magnitude() << '\t' << aoa_deg.Magnitude()  << '\t'
      << m.Magnitude()       << '\t' << ToString(m_mode)     << '\t'
      << tkg.Magnitude()     << '\t' << burnRate.Magnitude() << '\t'
      << Q.Magnitude()       << '\t' << M                    << '\t'
      << longG               << '\t' << pa.Magnitude()       << '\t'
      << a_dt.Magnitude()    << std::endl;
  }
  return cont;
}

//===========================================================================//
// Angle-of-Attack and Thrust Vector Elevation:                              //
//===========================================================================//
//        AoA   theta
std::pair<Angle,Angle> Ascent2::AoA(Time a_t, Angle a_psi) const
{
  // XXX: We do not check for sin(a_psi) < 0 here (ie switching from Ascent to
  // Descent motion),  as this would require a wider context than available in
  // this method. Such checks are performed in "ODECB".
  //
  // Also, the "retrograde" movement (|psi| > Pi/2, and thus omega < 0,  since
  // omega ~ Vhor ~ cos(psi)) may in general occur due to rounding errors or
  // integration instabilities. It is difficult to say precisely which values
  // of "psi" are accepatble and which are not, so we just log such events:
  //
  if (Abs(a_psi) > PI_2 && Base::m_os != nullptr && Base::m_logLevel >= 1)
#   pragma omp critical(Output)
    *Base::m_os << "Ascent2::AoA: WARNING: psi="
                << To_Angle_deg(a_psi).Magnitude() << " deg" << std::endl;

  Angle  aoa = 0.0_rad;

  switch (m_mode)
  {
  case FlightMode::Burn1:
  {
    // Arg: "tau" which is the time since Stage1 ignition (where the latter is
    // calculated via "m_T1", because the actual event-based "m_ignTime1"   is
    // not reached yet), so AoA1(tau=0)=0, ie @ launch;   0 <= tau <= T1:
    Time tIgn1 = m_cutOffTime1  - Base::m_T1;
    Time tau   = std::min(std::max(a_t - tIgn1,    0.0_sec), Base::m_T1);
    aoa        = std::max(tau * (m_aAoA1 * tau + m_bAoA1), 0.0_rad);
    if (aoa > 1.01 * m_maxAoA1)
    {
      OutputCtls();
      throw std::logic_error
            ("Burn1: MaxAoA1=" + std::to_string(m_maxAoA1.Magnitude()) +
             " exceeded: AoA=" + std::to_string(aoa      .Magnitude()));
    }
    break;
  }

  case FlightMode::Gap:
    break;

  case FlightMode::Burn2:
  {
    // Arg: "t" is time from Stage2 cut-off (same as orbital insertion time),
    // so it is just the main time "a_t": -T2 <= t <= 0; enforce those constr-
    // aints. NB: Formulas similar to those for Stage2, but with the INVERTED
    // sign of "b"; AoA2(t=0)=0, ie @ orbital insertion:
    assert(!IsPos(a_t));
    Time t = std::max(a_t, -m_T2);
    aoa    = std::max(t *  (m_aAoA2 * t - m_bAoA2), 0.0_rad);
    if (aoa > 1.01 * m_maxAoA2)
    {
      OutputCtls();
      throw std::logic_error
            ("Burn2: MaxAoA2=" + std::to_string(m_maxAoA2.Magnitude()) +
             " exceeded: AoA=" + std::to_string(aoa      .Magnitude()));
    }
    break;
  }
  default: ;
  }

  // XXX: In this case, we assume there is no chamber/nozzle gimbaling,  and
  // therefore, theta = psi + AoA; but "theta" is limited by Pi/2, otherwise
  // the Thrust vector will point towards negative polar angles. Yet theta < 0
  // is allowed under some circumstances (because Vr < 0 may be allowed):
  Angle theta = a_psi + aoa;

  if (theta > PI_2)
  {
    aoa   =  PI_2 - a_psi;
    theta =  PI_2;
  }
  else
  if (theta < -PI_2)
  {
    aoa   = -PI_2 - a_psi;
    theta = -PI_2;
  }
  return std::make_pair(aoa, theta);
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
Ascent2::AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const
{
  auto atm  = EAM::GetAtmConds(std::max(a_r - R, 0.0_km));
  double M  = Base::Mach(atm, a_v);
  assert(IsFinite(M) &&  M >= 0.0);

  // The Aerodynamic Force Main Term:
  ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * Base::m_crosS;

  // The Drag and Lift Coeffs, using the default model:
  double cD = LVAeroDyn::cD(M, a_AoA);
  double cL = LVAeroDyn::cL(M, a_AoA);

  return std::make_tuple(atm, cD * F, cL * F);
}

//===========================================================================//
// Propellant Burn Rate: May be variable:                                    //
//===========================================================================//
MassRate Ascent2::PropBurnRate
(
  Mass UNUSED_PARAM(a_m),
  LenK UNUSED_PARAM(a_h),
  Time a_t
)
const
{
  switch (m_mode)
  {
  case FlightMode::Burn1:
  {
    // Arg: "tau" is the time since Stage1 ignition  (so 0 <= tau <= T1); but
    // the latter is not yet reached, so calculate it using "T1":
    Time     tIgn1 = m_cutOffTime1 - Base::m_T1;
    Time     tau   = std::min(std::max(a_t  - tIgn1,    0.0_sec), Base::m_T1);
    MassRate res   =
      std::max(Base::m_burnRateI1  + (m_bMu1 + m_aMu1 * tau) * tau,
               MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));
    return res;
  }
  case FlightMode::Gap:
    return MassRate(0.0);

  case FlightMode::Burn2:
  {
    // Arg: "tau" is the time since Stage2 ignition  (so 0 <= tau <= T2); but
    // the latter is not yet reached, so calculate it using "T2":
    Time     tIgn2 = - m_T2;
    Time     tau   = std::min(std::max(a_t  - tIgn2,    0.0_sec), m_T2);
    MassRate res   =
      std::max(m_burnRateI2        + (m_bMu2 + m_aMu2 * tau) * tau,
               MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));
    return res;
  }
  default:
    return MassRate(0.0);
  }
}

//===========================================================================//
// Thrust:                                                                   //
//===========================================================================//
ForceK Ascent2::Thrust(MassRate a_burn_rate, Pressure a_p) const
{
  assert(!(IsNeg(a_burn_rate) || IsNeg(a_p)));
  switch (m_mode)
  {
    case FlightMode::Burn1:
    {
      assert(IsPos(a_burn_rate));
      // Here we take the static air pressure effects into account:
      // Isp for this Pressure:
      double   x   = double(a_p / EAM::P0);
      assert(0.0 <= x && x <= 1.0);
      Time     Isp = x   * Base::m_IspSL1 + (1.0 - x) * Base::m_IspVac1;
      return a_burn_rate * Isp  * g0K;
    }

    case FlightMode::Gap:
      assert(IsZero(a_burn_rate));
      return ForceK(0.0);

    case FlightMode::Burn2:
      // We assume the Vacuum mode for Stage2:
      assert(IsPos(a_burn_rate));
      return a_burn_rate * m_IspVac2 * g0K;

    default:
      return ForceK(0.0);
  }
}

//===========================================================================//
// Current LV Mass (incl the PayLoad and Fairing if attached):               //
//===========================================================================//
Mass Ascent2::LVMass(Base::StateV const& a_s, Time a_t) const
{
  // The mass of Propellants spent between "a_t" < 0 and the Orbital Insertion
  // instant (t=0):
  Mass spentPropMass = std::get<3>(a_s);
  assert(!IsNeg(spentPropMass));

  // Initialise "m" to the mass @ Orbital Insertion (incl the unspent Stage2
  // Propellant):
  Mass m = m_emptyMass2 + m_unSpendable2 + m_payLoadMass;

  // In any case, ADD the SpentPropMass (between "a_t" and Orbital Insertion):
  m += spentPropMass;

  // If Stage1 has NOT separated yet, add its Empty and PropRemnants:
  if (m_mode != FlightMode::Burn2 && m_mode != FlightMode::Gap)
    m += Base::m_emptyMass1 + Base::m_unSpendable1;

  // Add the FairingMass if the Fairing has NOT separated @  "a_t". In the
  // Bwd time, it means that "m_fairingSepTime" is known and "a_t" is below
  // it:
  if (IsFinite(m_fairingSepTime) && a_t < m_fairingSepTime)
    m += m_fairingMass;

  // All Done:
  assert(IsPos(m));
  return m;
}

//===========================================================================//
// "OutputCtls":                                                             //
//===========================================================================//
void Ascent2::OutputCtls() const
{
  if (Base::m_os == nullptr || Base::m_logLevel < 1)
    return;

# pragma omp critical(Output)
  {
    // Here "t" is the time from Stage2 cut-off (so t=-T2..0):
    *Base::m_os << "# T2   := " << m_T2.Magnitude() << ';'    << std::endl;
    *Base::m_os << "# 0  <= tau <= T2"                        << std::endl;
    *Base::m_os << "#-T2 <= t   <= 0"                         << std::endl;
    *Base::m_os << "# AoA2 := t * ("   << m_aAoA2.Magnitude() << " * t - ("
                << m_bAoA2.Magnitude() << ")); "              << std::endl;
    *Base::m_os << "# mu2  := " << m_burnRateI2.Magnitude()   <<     " + ("
                << m_bMu2.Magnitude()  << ") * tau + ("
                << m_aMu2.Magnitude()  << ") * tau^2;"        << std::endl;

    // Here "t" is the time from Stage1 cut-off (so t=-T1..0), and "tau" is
    // the time since Stage1 ignition (NB: the (-b) coeff at "tau" for AoA),
    // so tau=0..T1:
    *Base::m_os << "# T1   := " << Base::m_T1.Magnitude() << ';'  << std::endl;
    *Base::m_os << "# 0  <= tau <= T1"                        << std::endl;
    *Base::m_os << "# AoA1 := tau * (" << m_aAoA1.Magnitude() << " * tau + ("
                << m_bAoA1.Magnitude() << ")); "              << std::endl;
    *Base::m_os << "# mu1  := " << Base::m_burnRateI1.Magnitude() << " + ("
                << m_bMu1.Magnitude()  << ") * tau + ("
                << m_aMu1.Magnitude()  << ") * tau^2;"        << std::endl;
  }
}
}
// End namespace SpaceBallistics
