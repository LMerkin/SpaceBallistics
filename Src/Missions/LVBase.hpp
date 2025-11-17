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
  Angle  psi    = Angle(ATan2(a_Vr, a_Vhor));
  double cosPsi = double(a_Vhor / V);
  double sinPsi = double(a_Vr   / V);

  // "AoA"  : the Angle-of-Attack;
  // "theta": elevation angle of the Thrust vector
  // ("psi" induces constraints on both).
  // If there is no gimbaling of thrust chambers / nozzles, then the Thrust
  // vector ctl  is only achieved via the AoA, and thus  theta = psi + AoA:
  //
  auto [aoa, theta] = ToDer()->AoA(a_t, psi);

  // AeroDynamic Drag and Lift:
  auto [atm, drag, lift]  = ToDer()->AeroDynForces(a_r, V, aoa);
  assert(!(IsNeg(drag) || IsNeg(lift)));

  // Propellant Burn Rate (>= 0) and Thrust:
  MassRate burnRate = ToDer()->PropBurnRate(a_m, a_r - R, a_t);
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
  // Components of "u" in the (radius-vector, normal-to-radius-vector) frame:
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
LVBase<Derived>::ODERHS   (StateV const& a_s, Time a_t, bool a_is_ascent) const
{
  assert((a_is_ascent && !IsPos(a_t)) || (!a_is_ascent && !IsNeg(a_t)));

  // NB:
  // (*) In the RHS evaluation, r <= R is allowed, as it does not cause any
  //     singularities by itself; but we detect this condition in the Call-Back,
  //     which means that the integration is over (successfully or otherwise);
  // (*) "spentPropMass" may sometimes become < 0 due to rounding errors, so
  //     guard against that:
  LenK   r                = std::get<0>(a_s);
  VelK   Vr               = std::get<1>(a_s);
  AngVel omega            = std::get<2>(a_s);
  Mass   spentPropMass    = std::max(std::get<3>(a_s), 0.0_kg);
  Angle  phi              = std::get<4>(a_s);
  assert(IsPos(r));
  // The "horizontal" velocity (orthogonal to the radius-vector):
  VelK   Vhor             = r * omega / 1.0_rad;

  // Singular Point detection criteria:  NB: near the Singularity,  "Vhor"
  // approaches 0 much faster than "Vr". If both |Vhor| and |Vr| are below
  // their resp thresholds, we assume that we are near the Singular Point.
  // In that case, we disregard "Vhor" and use "Vr" only.
  // XXX: STRICTLY SPEAKING, V=0 does not always imply a singularity. It may
  // theoretically happen, for example, that the Thrust vector is defined in-
  // dependently of "psi" (which requires V > 0 to be determined),  and then
  // the ODE RHS is well-defined even @ V=0. However, in practice V=0 cannot
  // occur under any normal flight conditions, so it is a singularity indeed:
  //
  if (Abs(Vhor) < SingVhor && Abs(Vr) < SingV)
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

  // NB: the "spentPropMass" derivative depends on the integration mode:
  // (Ascent = Bwd, Descent = Fwd):
  //
  // In the Descent(Fwd) mode, spentPropMass = int_0^t burnRate(t') dt',
  // therefore
  //   dot(spentPropMass)  =  burnRate  >= 0;
  // in the Ascent (Bwd) mode, spentPropMass = int_t^0 burnRate(t') dt'
  // (where t <= 0), and thus
  //   dot(spentPropMass)  = -burnRate  <= 0,
  // but in both cases, spentPropMass >= 0:
  //
  MassRate spmDot = a_is_ascent ? -burnRate : burnRate;

  assert(( a_is_ascent && !IsPos(spmDot)) ||
         (!a_is_ascent && !IsNeg(spmDot)));

  // The result:
  return std::make_tuple(Vr, r2Dot, omegaDot, spmDot, omega);
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
LVBase<Derived>::LocateSingularPoint
(
  NearSingularityExn const& a_nse,
  bool                      a_is_ascent
)
const
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  // TOLERANCES:
  constexpr LenK HTol = 0.05_km;
  constexpr VelK VTol = VelK(0.01);

  // If we got here, the horizontal velocity is considered to be negligible;
  // only use the radial one which is equal to the total velocity:
  LenK  r1  = a_nse.m_r;
  LenK  h1  = r1 - R;
  VelK  Vr1 = a_nse.m_Vr;
  Time  t1  = a_nse.m_t;
  LenK  LS  = R * double(a_nse.m_phi);

  // The Curr LV Mass (via a synthetic "s1"):
  StateV s1    { r1, Vr1, AngVel(0.0), a_nse.m_spentPropMass, a_nse.m_phi };
  Mass   m1  = ToDer()->LVMass(s1, t1);
  assert(IsPos(m1));

  // Constraints:
  Pressure maxQ     = ToDer()->m_maxQ;
  Pressure sepQ     = ToDer()->m_sepQ;
  double   maxLongG = ToDer()->m_maxLongG;

  // NB: "Vr1"  is assumed to be small; typically:
  // (*) for the Ascent-To-Orbit (Bwd integration), Vr1 > 0;
  // (*) for the Descent (RTLS)  (Fwd integration), Vr1 < 0:
  // XXX: First, if "Vr1" is really small, we simply reset it to 0:
  if (Abs(Vr1) < VTol)
    Vr1 = VelK(0.0);

  // Otherwise: Check for Vr < 0:
  // In the Ascent mode,  Vr < 0 may occur (eg near the Orbital Insertion), but
  // definitely not near the Singularity;   in the Descent mode, Vr < 0 is nor-
  // mal:
  if (a_is_ascent && IsNeg(Vr1))
  {
    if (m_os != nullptr)
#     pragma omp critical(Output)
      *m_os << "# LVBase::LocateSingularPoint: ERROR: Ascent, but Vr1="
            << Vr1.Magnitude() << " km/sec" << std::endl;
    // Return Error:
    return RunRes();
  }
  // Thus: at this point, Ascent => Vr > 0:
  assert(!a_is_ascent || IsPos(Vr1));

  // Check for a degenerate "h1": it should not happen, but may:
  if (!IsPos(h1))
  {
    if (m_os != nullptr && m_logLevel >= 1 && h1 < -HTol)
#     pragma omp critical(Output)
      *m_os << "# LVBase::LocateSingularPoint: WARNING: h1=" << h1.Magnitude()
            << " km; reset to 0" << std::endl;

    // Assuming that the negative vals of "h1" if occur, are anyway small, so
    // reset them to 0:
    r1  = R;
    h1  = 0.0_km;
  }
  assert(r1 >= R && !IsNeg(h1));

  // The "equivalent" velocity @ h=0, using the Energy Integral:
  VelK V0 = SqRt(Sqr(Vr1) + K * (2.0 / R - 2.0 / r1));

  //-------------------------------------------------------------------------//
  // Forces and Acceleration:                                                //
  //-------------------------------------------------------------------------//
  // Assume that only the Gravitational Force and the Thrust are applicable
  // (the AeroDynamic Force can be neglected because the velocity is near-0,
  // yet we still need the "atm" to adjust the Thrust),  and the LV Mass is
  // CONSTANT during this final short interval. The assumption LVMass=const
  // is valid up to O(Vr1), as can be seen from the exact solution (using the
  // LambertW function):
  //
  AccK  g1 = K / Sqr(r1);
  auto atm = EAM::GetAtmConds(h1);

  // Propellant BurnRate (>= 0) and Thrust:
  MassRate burnRate1 = ToDer()->PropBurnRate(m1, h1, t1);
  assert(!IsNeg(burnRate1));

  // We will assume that the Thrust vector is pointing strictly UpWards
  // (psi = theta = Pi/2, alpha = 0):
  ForceK   thrust1   = ToDer()->Thrust(burnRate1, std::get<0>(atm));
  assert(!IsNeg(thrust1) && IsZero(thrust1) == IsZero(burnRate1));

  // The effective (radial) acceleration, estimated "manually" (ie NOT re-using
  // the "ODERHS" method):
  // (*) Since near the Singularity omega=0, there are no "kinematic" terms in
  //     the acceleration;
  // (*) and since V=~0, we drop the aerodynamic terms;
  // (*) XXX: Furthermore, we ASSUME that the Thrust vector points in the radi-
  //     us-vector direction, ie theta=psi=Pi/2 here, AoA=0 but is irrelevant,
  //     otherwise we would not be able to use the following simple 1D computa-
  //     tions to locate the Singular Point!..
  //
  // Then the radial acceleration is:
  AccK acc1 = thrust1 / m1 - g1;

  //-------------------------------------------------------------------------//
  // Edge Cases:                                                             //
  //-------------------------------------------------------------------------//
  if (IsZero(h1))
    return RunRes(RunRC::ZeroH,       t1, LS, Abs(Vr1), Abs(acc1), m1,
                  maxQ, sepQ, maxLongG);

  // If Vr1=0, then we are at the Singular Point -- but with h1 > 0, we must
  // return the equivalent "V0":
  if (IsZero(Vr1))
    return RunRes(RunRC::Singularity, t1, LS, V0,  Abs(acc1), m1,
                  maxQ, sepQ, maxLongG);

  if (!IsPos(acc1))
  {
    // Then the LV is Falling Back to Earth, and the Singular Point (where V=0)
    // is not reachable:
    if (a_is_ascent)
    {
      // During the Ascent, this is definitely an error cond, which is very un-
      // likely to occur:
      if (m_os != nullptr)
#       pragma omp critical(Output)
        *m_os << "# LVBase::LocateSingularPoint: ERROR: UnReachable: Acc="
              << double(acc1 / g1) << " g, Vr="      << Vr1.Magnitude()
              << " km/sec, h="     << h1.Magnitude() << " km" << std::endl;
      // Return Error:
      return RunRes();
    }
    else
    {
      // During the Descent, this may happen if we got a too small thrust; this
      // simply means that there will be hard landing, so  return a "synthetic"
      // "ZeroH" condition:
      if (m_os != nullptr && m_logLevel >= 1)
#       pragma omp critical(Output)
        *m_os << "# LVBase::LocateSingularPoint: WARNING: Descending, and Acc="
              << double(acc1 / g1) << " g, Vr="      << Vr1.Magnitude()
              << " km/sec, h="     << h1.Magnitude() << " km" << std::endl;

      // Estimate the Fall Time (w/o the aerodynamic forces, and assuming the
      // constant Thrust but still 

      // (*) "V0" is approximated by the Energy Integral (see above);
      // (*) the final LVMass is taken to be the curr one -- this is OK in the
      //     Descent mode;
      // (*) XXX: We do not compute the Fall Time, so use NAN (not "t1"):
      //
      return RunRes(RunRC::ZeroH, Time(NAN), LS, V0, Abs(acc1), m1,
                    maxQ,   sepQ, maxLongG);
    }
  }
  // So: The remaining Generic Case:
  assert(r1 > R && IsPos(h1) && IsPos(acc1));

  //-------------------------------------------------------------------------//
  // Remaining Time and Distance to the Singular Point:                      //
  //-------------------------------------------------------------------------//
  Time tau = -Vr1            / acc1;    // Of any sign!
  LenK dr  = -0.5 * Sqr(Vr1) / acc1;    // Consistent with the exact solution!
  assert(!(IsPos(dr)  || IsZero(tau)));

  // If follows from the above that for Ascent, we always have Vr1 > 0 and
  // acc1 > 0, so tau < 0, which is perfectly correct in the Bwd Integration
  // mode:   IsAscent => tau < 0:
  assert(!a_is_ascent || IsNeg(tau));

  // In the Descent mode (Fwd Integration), we should normally have tau > 0.
  // The opposite case requires a special handling:
  if (!a_is_ascent && IsNeg(tau))
  {
    // Since acc1 > 0 and tau < 0, then Vr1 > 0 in the Descent mode, which is
    // wrong (somehow, we have started moving away up again):
    assert(IsPos(Vr1));

    if (m_os != nullptr && m_logLevel >= 1)
#     pragma omp critical(Output)
      *m_os << "# LVBase::LocateSingularPoint: WARNING: Descending, but Vr1="
            << Vr1.Magnitude() << " km/sec, h1=" << std::endl;

    // XXX: Still, since acc1 < 0, we assume that we will eventually fall back
    // to Earth with the velocity "V0", although the fall time is UNKNOWN (NAN):
    //
    return RunRes(RunRC::ZeroH, Time(NAN), LS, V0, Abs(acc1), m1,
                  maxQ,  sepQ,  maxLongG);
  }
  // Thus:
  assert(a_is_ascent == IsNeg(tau));

  //-------------------------------------------------------------------------//
  // Finally: State and Time of the Singular Point:                          //
  //-------------------------------------------------------------------------//
  Time tS  = t1 + tau;
  LenK rS  = r1 + dr;
  LenK hS  = rS - R;

  // We should have hS >= 0; negative vals will be converted to 0:
  if (IsNeg(hS))
  {
    if (m_os != nullptr && m_logLevel >= 1 && hS < -HTol)
#     pragma omp critical(Output)
      *m_os << "# LVBase::LocateSingularPoint: WARNING: Got hS="
            << hS.Magnitude() << " km; reset to 0" << std::endl;
    rS = R;
    hS = 0.0_km;
  }
  // If we got here: Singular Point has been located successfully:
  // NB: "propS" is the total SpentPropMass at the Singular Point;
  // the following formula is valid for both Bwd and Bwd modes, since tau >= 0:
  Mass propS = a_nse.m_spentPropMass + burnRate1 * Abs(tau);

  StateV singS {rS, VelK(0.0), AngVel(0.0), propS, a_nse.m_phi};
  Mass mS    = ToDer()->LVMass(singS, tS);

  // Another way of computing "mS":
  Mass mSalt = m1 + burnRate1 * Abs(tau);
  if (!mS.ApproxEquals(mSalt) && m_os != nullptr)
#   pragma omp critical(Output)
    *m_os << " LVBase::LocateSingularPoint: ERROR: Inconsistency: mS="
          << mS.Magnitude() << ", mSalt=" << mSalt.Magnitude() << std::endl;

  if (m_os != nullptr && m_logLevel >= 3)
  {
#   pragma omp critical(Output)
    *m_os << "# SingularPoint Located: t1="     <<   t1.Magnitude()
          << " sec, tau="   << tau.Magnitude()  << " sec, tS="
          << tS.Magnitude() << " sec, m1="      <<   m1.Magnitude()
          << " kg, mS="     << mS.Magnitude()   << " kg, hS="
          << (rS - R).Magnitude()               << " km, LS="
          << LS      .Magnitude()               << " km, "
          << ToDer()->ToString(ToDer()->m_mode) << std::endl;
  }
  // Convert "rS" into the velocity @ H=0 (to allow for uniform treatment of
  // constraints) using the Energy Integral:  Adjust "V0" using "rS":
  assert(rS >= R);
  V0 = SqRt(K * (2.0 / R - 2.0 / rS));

  // Re-calculate (approximately) the final (radial) acceleration:
  AccK accS = thrust1 / mS - g1;

  // The final result:
  return RunRes(RunRC::Singularity, tS, LS, V0, Abs(accS), mS,
                maxQ, sepQ, maxLongG);
}

//===========================================================================//
// "PostProcessRun":                                                         //
//===========================================================================//
template<typename Derived>
LVBase<Derived>::RunRes
LVBase<Derived>::PostProcessRun(StateV const& a_sT, Time a_T, bool a_is_ascent)
const
{
  //-------------------------------------------------------------------------//
  // The Final State:                                                        //
  //-------------------------------------------------------------------------//
  LenK     rEnd      = std::get<0>(a_sT);
  LenK     hEnd      = rEnd - R;
  VelK     VrEnd     = std::get<1>(a_sT);
  AngVel   omegaEnd  = std::get<2>(a_sT);
  Angle    phiEnd    = std::get<4>(a_sT);
  VelK     VhorEnd   = rEnd * omegaEnd / 1.0_rad;
  auto     VEnd2     = Sqr(VrEnd) + Sqr(VhorEnd);
  Mass     mEnd      = ToDer()->LVMass(a_sT, a_T);
  LenK     lEnd      = R * double(phiEnd);

  Pressure maxQ      = ToDer()->m_maxQ;
  Pressure sepQ      = ToDer()->m_sepQ;
  double   maxLongG  = ToDer()->m_maxLongG;

  // Ensure we are above the Earth surface:
  if (rEnd < R)
  {
    if (hEnd < -0.5_km && m_os != nullptr && m_logLevel >= 1)
#     pragma omp critical(Output)
      *m_os  << "# LVBase::PostProcessRun: WARNING: Arrived at h="
             << hEnd.Magnitude() << " km" << std::endl;
    rEnd = R;
    hEnd = 0.0_km;
  }

  //-------------------------------------------------------------------------//
  // Final Acceleration:                                                     //
  //-------------------------------------------------------------------------//
  // XXX: Unlike "LocateSingularPoint", here   we CANNOT in general assume that
  // omega=0, omegaDot=0 and the Thrust vector is pointing in the radius-vector
  // direction; so we have to invoke the "ODERHS" and compute the "accEnd" from
  // it:
  DStateV  dsT         = ODERHS(a_sT, a_T, a_is_ascent);
  AccK     r2DotEnd    = std::get<1>(dsT);
  AngAcc   omegaDotEnd = std::get<2>(dsT);

  // Radial and Horizontal Acceleration:
  AccK     accRadEnd   = r2DotEnd - rEnd * Sqr(omegaEnd / 1.0_rad);
  AccK     accHorEnd   = (rEnd * omegaDotEnd + 2 * VrEnd * omegaEnd) / 1.0_rad;

  // Total Final Acceleration:
  AccK     accEnd      = SqRt(Sqr(accRadEnd) + Sqr(accHorEnd));

  //-------------------------------------------------------------------------//
  // The Result:                                                             //
  //-------------------------------------------------------------------------//
  RunRC  rc;
  if (ToDer()->m_mode == Derived::FlightMode::UNDEFINED)
    // We have run out of propellant:
    rc = RunRC::FlameOut;
  else
  {
    // Then we should get hEnd==0 (because the only other possibility is
    // that the integration ran until tMin, which is extremely unlikely
    // (XXX: TODO: still handle this condition properly!):
    //
    if (hEnd > 0.5_km && m_os != nullptr && m_logLevel >= 1)
#     pragma omp critical(Output)
      *m_os  << "# LVBase::PostProcessRun: WARNING: Expected h=0 but got h="
             << hEnd.Magnitude() << " km" << std::endl;

    // Still, formally it is a "ZeroH" outcome (this is OK, as we still correct
    // the final velocity for "rEnd"):
    rc = RunRC::ZeroH;
  }
  assert(!IsNeg(hEnd) && rEnd >= R);

  // XXX: IMPORTANT: BEWARE: Convert the (rEnd, VEnd) into the equivalent Velo-
  // city @ h=0 using the Energy Integral, but use the actual "accEnd"! This is
  // OK for the moment:
  VelK   V0   = SqRt(VEnd2 + K * (2.0 / R - 2.0 / rEnd));

  return RunRes(rc, a_T, lEnd, V0, accEnd, mEnd, maxQ, sepQ, maxLongG);
}
}
// End namespace SpaceBallistics
