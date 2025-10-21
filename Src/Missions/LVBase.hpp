// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Missions/LVBase.hpp":                  //
//===========================================================================//
#include "SpaceBallistics/Missions/LVBase.h"

namespace SpaceBallistics
{
//===========================================================================//
// Non-Default Ctor:                                                         //
//===========================================================================//
template<typename Derived>
LVBase<Derived>::LVBase
(
  // Stage1:
  Mass            a_full_mass1,
  double          a_K1,               // PropMass1 / FullMass1
  double          a_prop_rem1,
  Time            a_Isp_sl1,
  Time            a_Isp_vac1, 
  ForceK          a_thrust_vac1,
  double          a_min_thrtl1,
  Len             a_diam,

  // Integration and Output Params:
  Time            a_ode_integr_step,
  std::ostream*   a_os,
  int             a_log_level
)
: // Const Stage1 Params:
  m_crosS         (0.25 * Pi<double> * Sqr(To_Len_km(a_diam))),
  m_K1            (a_K1),
  m_propRem1      (a_prop_rem1),
  m_IspSL1        (a_Isp_sl1),
  m_IspVac1       (a_Isp_vac1),
  m_minThrtL1     (a_min_thrtl1),

  // Stage1 Params which may be subject to Optimisation:
  m_fullMass1     (a_full_mass1),
  m_emptyMass1    (m_fullMass1 * (1.0 - m_K1)),
  m_propMass1     (m_fullMass1 * m_K1),
  m_unSpendable1  (m_propMass1 * m_propRem1),
  m_spendable1    (m_propMass1 - m_unSpendable1),
  m_thrustVacI1   (a_thrust_vac1),
  m_burnRateI1    (m_thrustVacI1 / (m_IspVac1 * g0K)),
  m_T1            (m_spendable1  / m_burnRateI1),

  // Integration and Output Params:
  m_odeIntegrStep (a_ode_integr_step),
  m_os            (a_os),             // May be NULL
  m_logLevel      (a_log_level)
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  if (!(0.0 < m_K1            && m_K1 < 1.0            &&
        IsPos(m_fullMass1)    && IsPos(m_emptyMass1)   &&
        IsPos(m_propMass1)    && IsPos(m_unSpendable1) &&
        IsPos(m_spendable1)   &&
        IsPos(m_IspVac1)      && IsPos(m_IspSL1)       &&
        m_IspSL1 < m_IspVac1  && IsPos(m_burnRateI1)   &&
        IsPos(m_thrustVacI1)  && IsPos(m_T1)           &&
        0.0 <= m_minThrtL1    && m_minThrtL1  <= 1.0   &&
        IsPos(m_crosS)        && IsPos(m_odeIntegrStep)))
    throw std::invalid_argument("LVBase::Ctor: Invalid LV/ODE Param(s)");
}

//===========================================================================//
// The Mach Number:                                                          //
//===========================================================================//
template<typename Derived>
double LVBase<Derived>::Mach(EAM::AtmConds const& a_atm, VelK a_v)
{
  assert(!IsNeg(a_v));

  // The Speed of Sound:
  Vel A = std::get<3>(a_atm);
  return
    IsZero(A)
    ? Inf<double>                 // We are above the atmosphere
    : double(To_Len_m(a_v) / A);  // Genertic case
}

//===========================================================================//
// "NonGravForces":                                                          //
//===========================================================================//
// Common Helper for "ODERHS" and "ODECB":
//
template<typename Derived>
void LVBase<Derived>::NonGravForces
(
  Time           a_t,
  LenK           a_r,
  VelK           a_Vr,
  VelK           a_Vhor,
  Mass           a_m,
  VelK*          a_V,
  Angle*         a_psi,
  Angle*         a_aoa,
  EAM::AtmConds* a_atm,
  MassRate*      a_burn_rate,
  ForceK*        a_thrust,
  double         a_lv_axis[2],  // In the (r, normal-to-r) frame
  AccK           a_ng_acc [2],  // ditto
  double*        a_long_g       // Longitudinal G
)
const
{
  assert(IsPos(a_r)          && IsPos(a_m));
  assert(a_V      != nullptr && a_psi       != nullptr && a_aoa    != nullptr &&
         a_atm    != nullptr && a_burn_rate != nullptr && a_thrust != nullptr &&
         a_ng_acc != nullptr && a_lv_axis   != nullptr && a_long_g != nullptr);

  // The Absolute Velocity, should be bounded away from 0:
  VelK   V  = SqRt(Sqr(a_Vr) + Sqr(a_Vhor));
  assert(IsPos(V));

  // Trajectory Inclination Angle. NB:
  // We have a_Vhor >= 0, and normally a_Vr >= 0, so in the generic case,
  // 0  <= psi <=  Pi/2;
  // at start, psi=Pi/2; at orbital insertion, psi=0;
  // sinPsi may be <  0 if we Fall Back to Earth (in that case Psi = -Pi/2);
  // cosPsi may be <  0 if we are returning to the launch site;
  // XXX: Yet, depending on the mission profile, Vr < 0 (or equiv. sinPsi < 0)
  // may or may not be allowed; but this is decided in the Derived::ODERHS, not
  // here:
  Angle  psi     = Angle(ATan2(a_Vr, a_Vhor));
  double cosPsi  = double(a_Vhor / V);
  double sinPsi  = double(a_Vr   / V);

  // "AoA"  : the Angle-of-Attack;
  // "theta": elevation angle of the Thrust vector
  // ("psi" induces constraints on both).
  // If there is no gimbaling of thrust chambers / nozzles, then the Thrust vec-
  // tor ctl  is only achieved via the AoA, and thus theta = psi + AoA:
  //
  auto [aoa, theta] = ToDer()->AoA(a_t, psi);
  assert(!IsNeg(aoa) && aoa < PI_2 && !IsNeg(theta) && theta <= PI_2);

  // AeroDynamic Drag and Lift:
  auto [atm, drag, lift]  = ToDer()->AeroDynForces(a_r, V, aoa);
  assert(!(IsNeg(drag) || IsNeg(lift)));

  // Propellant Burn Rate (>= 0) and Thrust:
  MassRate burnRate = ToDer()->PropBurnRate(a_t);
  assert(!IsNeg(burnRate));

  Pressure p        = std::get<0> (atm);
  ForceK   thrust   = ToDer()->Thrust(burnRate, p);
  assert(!IsNeg(thrust) && IsZero(thrust) == IsZero(burnRate));

  // "v" unit vector: in the velocity direction; "psi" introduced above is the
  // velocity elevation over the horizon   (ie over the positive normal to the
  // radius-vector):
  double v[2]  { sinPsi,  cosPsi };
  //               r    normal-to-r
  // "u" unit vector: normal to the velocity (in the "up" direction):
  // Components of "u" in the (radius-vector, normal-to-radius-vector)
  // frame:
  double u[2]  { cosPsi, -sinPsi };
  //               r    normal-to-r

  // Thrust components (in the same frame):
  double cosTheta = Cos(theta);
  double sinTheta = Sin(theta);
  // If there is no gimbaling and theta = psi + AoA, then also
  // cos(theta) = v[1] * cosA + u[1] * sinA;
  // sin(theta) = v[0] * cosA + u[0] * sinA;
  //
  ForceK T[2]
  {
    sinTheta * thrust,   // r
    cosTheta * thrust    // normal-to-r
  };

  // AeroDynamic Drag (-v) and Lift (+u) components:
  //                   r        normal-to-r
  ForceK AD[2] { -v[0] * drag, -v[1] * drag };
  ForceK AL[2] {  u[0] * lift,  u[1] * lift };

  // Acceleration components due to the above forces:
  a_ng_acc [0] = (T[0] + AD[0] + AL[0]) / a_m;  // r
  a_ng_acc [1] = (T[1] + AD[1] + AL[1]) / a_m;  // normal-to-r

  // The components of the logitudinal (X) axis of the LV in the
  // (r, normal_to_r) frame:
  // BEWARE: here we really have to use (psi + AoA), which may be different
  // from "theta":
  if (psi + aoa == theta)
  {
    a_lv_axis[0] = sinTheta;
    a_lv_axis[1] = cosTheta;
  }
  else
  {
    a_lv_axis[0] = Sin(psi + aoa);
    a_lv_axis[1] = Cos(psi + aoa);
  }
  // The Longitudinal G:
  // It is a projection of "ngAcc" to the main LV axis (XXX: verify whether
  // this expr is also correct in case of Falling Back to Earth --  but we
  // only need the abs val). In the (r, normal_to_r) frame, the longitudinal
  // LV axis is given by the (psi + aoa) angle as above:
  *a_long_g    =
    double(Abs(a_lv_axis[0] * a_ng_acc[0] + a_lv_axis[1] * a_ng_acc[1]) / g0K);

  // Return other "physical" variables as well:
  *a_V         = V;
  *a_psi       = psi;
  *a_aoa       = aoa;
  *a_atm       = atm;
  *a_burn_rate = burnRate;
  *a_thrust    = thrust;
}

//===========================================================================//
// ODE RHS:                                                                  //
//===========================================================================//
template<typename Derived>
typename LVBase<Derived>::DStateV
LVBase<Derived>::ODERHS   (StateV const& a_s, Time a_t) const
{
  // NB: In the RHS evaluation, r <= R is allowed, as it does not cause any
  // singularities by itself; but we detect this condition in the Call-Back,
  // which means that the integration is over (successfully or otherwise):
  LenK   r                = std::get<0>(a_s);
  VelK   Vr               = std::get<1>(a_s);
  AngVel omega            = std::get<2>(a_s);
  Mass   spentPropMass    = std::get<3>(a_s);
  Angle  phi              = std::get<4>(a_s);
  assert(!IsPos(a_t) && IsPos(r));
  // The "horizontal" velocity (orthogonal to the radius-vector):
  VelK   Vhor             = r * omega / 1.0_rad;

  // Singular Point detection criteria:  NB: near the Singularity,  "Vhor"
  // approaches 0 much faster than "Vr". If both |Vhor| and |Vr| are below
  // their resp thresholds, we assume that we are near the Singular Point.
  // In that case, we disregard "Vhor" and use "Vr" only:
  //
  if (Abs(Vhor) < SingVhor && Abs(Vr) < SingVr)
    throw NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  // Generic Case:
  // The curr mass:
  Mass   m  = ToDer()->LVMass(a_s, a_t);
  assert(IsPos(m));

  // Compute the Non-Gravitational Acceleration Components:
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

  // The RHS components:
  AccK r2Dot =
      r * Sqr(omega / 1.0_rad)   // "Kinematic"       term
    - K / Sqr(r)                 // Gravitational     term
    + ngAcc[0];                  // Non-Gravitational accs

  AngAcc omegaDot =
    (
      - 2.0 * Vr * omega         // "Kinematic"       term
      + ngAcc[1] * 1.0_rad       // Non-Gravitational accs
    )
    / r;

  // The result: NB: the last derivative is <= 0, since the corresp "StateV"
  // components is int_t^0 burnRate(t') dt', ie "t" is the LOWER integration
  // limit:
  return std::make_tuple(Vr, r2Dot, omegaDot, -burnRate, omega);
}

//===========================================================================//
// "LocateSingularPoint":                                                    //
//===========================================================================//
// The final stage of integration, where we assume omega=0 (a purely vert-
// ical motion)   and integrate the simplified ODEs ANALYTICALLY to avoid
// numerical instabilities in the viciniy of the singular point.
// Returns (singH, singT)  if the singular point has been found, otherwise
// "nullopt":
//
template<typename Derived>
LVBase<Derived>::RunRes
LVBase<Derived>::LocateSingularPoint(NearSingularityExn const& a_nse) const
{
  // If we got here, the horizontal velocity is considered to be negligible;
  // only use the radial one which is equal to the total velocity:
  LenK  r1  = a_nse.m_r;
  LenK  h1  = r1 - R;
  VelK  Vr1 = a_nse.m_Vr;
  Time  t1  = a_nse.m_t;

  // Check for the following degenerate conditions: They should not happen,
  // but may:
  if (!(IsPos(h1) && IsPos(Vr1)) && m_os != nullptr && m_logLevel >= 1)
    *m_os << "Ascent2::LocateSingularPoint: WARNING: h="
          << h1.Magnitude() << " km, Vr=" << Vr1.Magnitude() << " km/sec"
          << std::endl;
  // Assuming that the negative vals of "h1" or "Vr1", if occur, are anyway
  // small, round them up:
  if (IsNeg(h1))
  {
    r1  = R;
    h1  = 0.0_km;
  }
  if (IsNeg(Vr1))
    Vr1 = VelK(0.0);

  // Assume that only the Gravitational Force and the Thrust are applicable
  // (the AeroDynamic Force can be neglected because the velocity is near-0,
  // yet we still need the air pressure to adjust the Thrust), and the LV
  // Mass is CONSTANT at this final short interval; yet we need the "atm" for
  // Thrust conputation:
  AccK  g1 = K / Sqr(r1);
  auto atm = EAM::GetAtmConds(h1);

  // Propellant BurnRate (>= 0) and Thrust:
  MassRate burnRate1 = ToDer()->PropBurnRate(t1);
  assert(!IsNeg(burnRate1));

  ForceK   thrust1   = ToDer()->Thrust(burnRate1, std::get<0>(atm));
  assert(!IsNeg(thrust1) && IsZero(thrust1) == IsZero(burnRate1));

  // The Curr LV Mass (via a synthetic "s1"):
  StateV s1    { r1, Vr1, AngVel(0.0), a_nse.m_spentPropMass, a_nse.m_phi };
  Mass   m1  = ToDer()->LVMass(s1, t1);
  assert(IsPos(m1));

  Pressure maxQ     = ToDer()->m_maxQ;
  Pressure sepQ     = ToDer()->m_sepQ;
  double   maxLongG = ToDer()->m_maxLongG;

  // The (constant) acceleration we will use:
  AccK   acc1 = thrust1 / m1 - g1;
  if (!IsPos(acc1))
  {
    // Then the LV is Falling Back to the pad, and the singular point is not
    // reachable; this is because we have arrived at the mass "m1"  which is
    // too large:
    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# Ascent2::LocateSingularPoint: UnReachable: Acc="
            << double(acc1 / g1) << " g" << std::endl;
    return RunRes
          {RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN),
           Pressure(NAN), Pressure(NAN), NAN};
  }

  // Otherwise: Remaining Bwd Time and Distance to the singular point:
  Time tau = Vr1            / acc1;
  LenK dr  = 0.5 * Sqr(Vr1) / acc1;
  assert(!(IsNeg(tau) || IsNeg(dr)));

  // Finally: State and Time of the singular point:
  LenK rS    = r1 - dr;
  LenK hS    = rS - R;

  // We should have hS >= 0; slightly negative vals will be rounded up to 0,
  // but larger negative vals are treated as errors:
  if (hS < -0.5_km)
  {
    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# Ascent2::LocateSingularPoint: Got hS=" << hS.Magnitude()
            << " km" << std::endl;
    return RunRes
          {RunRC::Error,  Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN),
           Pressure(NAN), Pressure(NAN),  NAN};
  }
  else
  if (IsNeg(hS))
  {
    rS = R;
    hS = 0.0_km;
  }

  // If we got here: Singular Point has been located successfully:
  Mass propS = a_nse.m_spentPropMass + burnRate1 * tau;
  StateV singS {rS, VelK(0.0), AngVel(0.0), propS, a_nse.m_phi};
  Time tS    = t1 - tau;
  Mass mS    = ToDer()->LVMass(singS, tS);
  LenK LS    = R * double(a_nse.m_phi);

  if (m_os != nullptr && m_logLevel >= 2)
  {
    *m_os << "# SingularPoint Located: t1="     <<     t1.Magnitude()
          << " sec, tau="   << tau.Magnitude()  << " sec, tS="
          << tS.Magnitude() << " sec, m1="      <<     m1.Magnitude()
          << " kg, mS="     << (m1 + burnRate1 * tau)    .Magnitude()
          << " kg, mSalt="  << mS.Magnitude()   << " kg, hS="
          << (rS - R).Magnitude()               << " km, LS="
          << LS      .Magnitude()               << " km, "
          << ToDer()->ToString(ToDer()->m_mode) << std::endl;
  }
  // Convert "rS" into the velocity @ H=0 (to allow for uniform treatment of
  // constraints) using the Energy Integral:
  assert(rS >= R);
  VelK V0 = SqRt(2.0 * K * (1.0 / R - 1.0 / rS));

  // The final result:
  return RunRes{RunRC::Singularity, tS, LS, V0, mS, maxQ, sepQ, maxLongG};
}

//===========================================================================//
// "PostProcessRun":                                                         //
//===========================================================================//
template<typename Derived>
LVBase<Derived>::RunRes
LVBase<Derived>::PostProcessRun(StateV const& a_sT, Time a_T) const
{
  // Check the final state and mode:
  //
  LenK     rEnd      = std::get<0>(a_sT);
  VelK     VrEnd     = std::get<1>(a_sT);
  AngVel   omegaEnd  = std::get<2>(a_sT);
  Angle    phiEnd    = std::get<4>(a_sT);
  auto     VEnd2     = Sqr(VrEnd) + Sqr(rEnd * omegaEnd / 1.0_rad);
  Mass     mEnd      = ToDer()->LVMass(a_sT, a_T);
  LenK     lEnd      = R * double(phiEnd);

  Pressure maxQ      = ToDer()->m_maxQ; 
  Pressure sepQ      = ToDer()->m_sepQ;
  double   maxLongG  = ToDer()->m_maxLongG;

  if (ToDer()->m_mode == Derived::FlightMode::UNDEFINED)
  {
    // We have run out of propellant while, presumably, still @ rEnd > R;
    // if we got rEnd < R, the diff should be very small:
    if (rEnd < R)
    {
      if (m_os != nullptr && m_logLevel >= 1)
        *m_os  << "# Ascent2::Run: Got mode=UNDEFINED but h="
               << (rEnd - R).Magnitude() << " km" << std::endl;
      rEnd = R;
    }
    assert(rEnd >= R);

    // Convert the (rEnd, VEnd) into the velocity:
    VelK V0 = SqRt(VEnd2 + 2.0 * K * (1.0 / R - 1.0 / rEnd));

    return RunRes{RunRC::FlameOut, a_T, lEnd, V0, mEnd, maxQ, sepQ, maxLongG};
  }
  else
  {
    // Then we should get hEnd==0 (because the only other possibility is
    // that the integration ran until tMin, which is extremely unlikely):
    LenK hEnd = rEnd - R;
    if  (hEnd > 0.5_km)
    {
      if (m_os != nullptr && m_logLevel >= 1)
        *m_os  << "# Ascent2::Run: Expected h=0 but got h="
               << hEnd.Magnitude() << " km" << std::endl;
      hEnd = 0.0_km; // Just reset it formally
    }
    return RunRes
          {RunRC::ZeroH, a_T, lEnd, SqRt(VEnd2), mEnd, maxQ, sepQ, maxLongG};
  }
}
}
// End namespace SpaceBallistics
