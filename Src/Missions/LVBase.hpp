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
  AccK           a_ng_acc [2]   // ditto
)
const
{
  assert(IsPos(a_r)          && !IsNeg(a_Vhor)         && IsPos(a_m));
  assert(a_V      != nullptr && a_psi       != nullptr && a_aoa    != nullptr &&
         a_atm    != nullptr && a_burn_rate != nullptr && a_thrust != nullptr &&
         a_ng_acc != nullptr && a_lv_axis   != nullptr);

  // The Absolute Velocity, should be bounded away from 0:
  VelK   V  = SqRt(Sqr(a_Vr) + Sqr(a_Vhor));
  assert(IsPos(V));

  // The Angle-of-Attack:
  Angle psi    = Angle(ATan2(a_Vr, a_Vhor)); // Trajectory Inclination Angle
  Angle aoa    = ToDer()->AoA(a_t, psi);
  assert(!IsNeg(aoa));
  double cosA  = Cos(aoa);
  double sinA  = Sin(aoa);

  // AeroDynamic Drag and Lift:
  auto [atm, drag, lift]  = ToDer()->AeroDynForces(a_r, V, aoa);
  assert(!(IsNeg(drag) || IsNeg(lift)));

  // Propellant Burn Rate (>= 0) and Thrust:
  MassRate burnRate = ToDer()->PropBurnRate(a_t);
  assert(!IsNeg(burnRate));

  Pressure p        = std::get<0> (atm);
  ForceK   thrust   = ToDer()->Thrust(burnRate, p);
  assert(!IsNeg(thrust) && IsZero(thrust) == IsZero(burnRate));

  // "v" unit vector: in the velocity direction: "Th" is the velocity elevation
  // of the hirizon (ie over the positive normal to the radius-vector):
  // v = [sinTh, cosTh];
  // at start, Th=Pi/2; at orbital insertion, Th=0:
  double sinTh = double(a_Vr   / V);
  double cosTh = double(a_Vhor / V);
  double v[2]  { sinTh,  cosTh };
  // NB:
  // sinTh may be < 0 if we Fall Back to Earth (in that case Th = -Pi/2); but
  // cosTh       >= 0 always (we do not move back in the polar angle):
  assert(cosTh   >= 0.0);

  // "u" unit vector: normal to the velocity (in the "up" direction):
  // Components of "u" in the (radius-vector, normal-to-radius-vector)
  // frame:
  double u[2]  { cosTh, -sinTh };

  // Thrust components (in the same frame):
  double sinApTh = v[0] * cosA + u[0] * sinA; // sin(A + Th)
  double cosApTh = v[1] * cosA + u[1] * sinA; // cos(A + Th) 
  ForceK T[2]
  {
    sinApTh * thrust,  // sin(A + Th) * thrust
    cosApTh * thrust   // cos(A + Th) * thrust
  };

  // AeroDynamic Drag (-v) and Lift (+u) components:
  ForceK AD[2] { -v[0] * drag, -v[1] * drag };
  ForceK AL[2] {  u[0] * lift,  u[1] * lift };

  // Acceleration components due to the above forces:
  a_ng_acc [0] = (T[0] + AD[0] + AL[0]) / a_m;
  a_ng_acc [1] = (T[1] + AD[1] + AL[1]) / a_m;

  // The components of the logitudinal (X) axis of the LV in the
  // (r, normal_to_r) frame:
  a_lv_axis[0] = sinApTh;
  a_lv_axis[1] = cosApTh;

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
LVBase<Derived>::ODERHS   (StateV const& a_s, Time a_t)
{
  // NB: In the RHS evaluation, r <= R is allowed, as it does not cause any
  // singularities by itself; but we detect this condition in the Call-Back,
  // which means that the integration is over (successfully or otherwise):
  LenK   r             = std::get<0>(a_s);
  VelK   Vr            = std::get<1>(a_s);
  AngVel omega         = std::get<2>(a_s);
  Mass   spentPropMass = std::get<3>(a_s);
  Angle  phi           = std::get<4>(a_s);
  assert(!IsPos(a_t) && IsPos(r));

  // The "horizontal" velocity (orthogonal to the radius-vector):
  VelK   Vhor          = r * omega / 1.0_rad;

  // NB: omega < 0 (or equivalently vHor < 0)  is qualitatively impossible,
  // because omega=0 implies omegaDot = 0; but it may occur due to a finite
  // integration step, so we have to control it manually:
  // If "Vhor" is not yet 0 but is below a certain positive threshold,  we
  // set omega=0 to avoid oscillations around the vertical:
  if (Vhor < SingVhor)
  {
    omega = AngVel(0.0);
    Vhor  = VelK  (0.0);
  }
  // If omega=0 and in addition "Vr" is below the threshold, we assume that
  // we are near the singular point:
  if (IsZero(omega) && Vr < SingVr)
    throw NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  // Generic Case:
  // Thus, omega < 0 cannot occur:
  assert(!IsNeg(omega));

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
  NonGravForces
    (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
     lvAxis, ngAcc);

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
std::optional<std::pair<typename LVBase<Derived>::StateV, Time>>
LVBase<Derived>::LocateSingularPoint(NearSingularityExn const& a_nse)
{
  // If we got here, the horizontal velocity is considered to be negligible;
  // only use the radial one:
  LenK  r1  = a_nse.m_r;
  LenK  h1  = r1 - R;
  VelK  Vr1 = a_nse.m_Vr;
  Time  t1  = a_nse.m_t;

  // Check for the following degenerate conditions: They should not happen,
  // but may:
  if (!(IsPos(h1) && IsPos(Vr1)) && m_os != nullptr && m_logLevel >= 2)
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

  // The (constant) acceleration we will use:
  AccK   acc1 = thrust1 / m1 - g1;
  if (!IsPos(acc1))
  {
    // Then the LV is Falling Back to the pad, and the singular point is not
    // reachable; this is because we have arrived at the mass "m1"  which is
    // too large:
    if (m_os != nullptr && m_logLevel >= 2)
      *m_os << "# Ascent2::LocateSingularPoint: UnReachable: Mass=Acc="
            << (double(acc1 / g1) - 1.0) << " g" << std::endl;
    return std::nullopt;
  }

  // Otherwise: Remaining Bwd Time and Distance to the singular point:
  Time tau = Vr1            / acc1;
  LenK dr  = 0.5 * Sqr(Vr1) / acc1;
  assert(!(IsNeg(tau) || IsNeg(dr)));

  // Finally: State and Time of the singular point:
  LenK rS    = r1 - dr;
  LenK hS    = rS - R;

  // We should have hS >= 0; slightly negative vals will be rounded up to 0:
  if (hS < -0.5_km)
  {
    if (m_os != nullptr)
      *m_os << "# Ascent2::LocateSingularPoint: Got hS=" << hS.Magnitude()
            << " km" << std::endl;
    return std::nullopt;
  }
  else
  if (IsNeg(hS))
  {
    rS = R;
    hS = 0.0_km;
  }

  Time tS    = t1 - tau;
  Mass propS = a_nse.m_spentPropMass + burnRate1 * tau;
  LenK LS    = R  * double(a_nse.m_phi);

  StateV singS {rS, VelK(0.0), AngVel(0.0), propS, a_nse.m_phi};

  if (m_os != nullptr && m_logLevel >= 2)
    *m_os << "# SingularPoint Located: t1="     <<     t1.Magnitude()
          << " sec, tau="   << tau.Magnitude()  << " sec, tS="
          << tS.Magnitude() << " sec, m1="      <<     m1.Magnitude()
          << " kg, mS="     << (m1 + burnRate1 * tau)    .Magnitude()
          << " kg, mSalt="  << ToDer()->LVMass(singS, tS).Magnitude()
          << " kg, hS="     << (rS - R).Magnitude()
          << " km, LS="     << LS      .Magnitude()
          << " km, "        << ToDer()->ToString(ToDer()->m_mode)
          << std::endl;

  return std::make_optional(std::make_pair(singS, tS));
}
}
// End namespace SpaceBallistics
