// vim:ts=2:et
//===========================================================================//
//                    "Src/Missions/Ascent2-Integr.cpp":                     //
//     Ascent-to-Orbit for a "Model" 2-Stage LV: Trajectory Integration      //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"
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

  // Logging Params:
  std::ostream*   a_os,
  int             a_log_level
)
: //-------------------------------------------------------------------------//
  // Over-All:                                                               //
  //-------------------------------------------------------------------------//
  // Const Over-All Params (not affected by Optimisation):
  m_maxStartMass  (a_max_start_mass),
  m_fairingMass   (a_fairing_mass),
  m_crosS         (0.25 * Pi<double> * Sqr(To_Len_km(a_diam))),

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
  m_maxAoA2       (To_Angle(a_max_aoa2)),

  // Non-Const Stage2 Params (may be updated in the course of Optimisation):
  m_fullMass2     ((m_maxStartMass  - m_fairingMass - m_payLoadMass) /
                   (1.0 + m_alpha1)),
  m_emptyMass2    (m_fullMass2   * (1.0 - m_K2)),
  m_propMass2     (m_fullMass2   * m_K2),
  m_unSpendable2  (m_propMass2   * m_propRem2),
  m_spendable2    (m_propMass2   - m_unSpendable2),
  m_thrustVacI2   (a_thrust_vac2),
  m_thrustMult2   (1.0),
  m_burnRateI2    (m_thrustVacI2 / (m_IspVac2 * g0K)),
  m_T2            (m_spendable2  / m_burnRateI2),

  //-------------------------------------------------------------------------//
  // Stage1:                                                                 //
  //-------------------------------------------------------------------------//
  // Const Stage1 Params (not affected by Optimisation):
  m_K1            (a_K1),
  m_propRem1      (a_prop_rem1),
  m_IspSL1        (a_Isp_sl1),
  m_IspVac1       (a_Isp_vac1),
  m_minThrtL1     (a_min_thrtl1),
  m_maxAoA1       (To_Angle(a_max_aoa1)),

  // Non-Const Stage1 Params (may be updated in the course of Optimisation):
  m_fullMass1     ((m_maxStartMass  - m_fairingMass - m_payLoadMass) /
                   (1.0 + m_alpha1) * m_alpha1),
  m_emptyMass1    (m_fullMass1      * (1.0 - m_K1)),
  m_propMass1     (m_fullMass1   * m_K1),
  m_unSpendable1  (m_propMass1   * m_propRem1),
  m_spendable1    (m_propMass1   - m_unSpendable1),
  m_thrustVacI1   (a_thrust_vac1),
  m_thrustMult1   (1.0),
  m_burnRateI1    (m_thrustVacI1 / (m_IspVac1 * g0K)),
  m_T1            (m_spendable1  / m_burnRateI1),

  //-------------------------------------------------------------------------//
  // Mission Params: Ininitiased in the Ctor Body:                           //
  //-------------------------------------------------------------------------//
  m_Rins          (0.0_km),
  m_Vins          (0.0),

  //-------------------------------------------------------------------------//
  // Flight Ctl Params:                                                      //
  //-------------------------------------------------------------------------//
  m_aMu2          (0.0),
  m_bMu2          (0.0),
  m_bHat2         (0.0),
  m_muHat2        (1.0),

  m_aAoA2         (0.0),
  m_bAoA2         (0.0),
  m_aAoAHat2      (0.0),
  m_bAoAHat2      (0.0),

  m_TGap          (0.0_sec),

  m_aMu1          (0.0),
  m_bMu1          (0.0),
  m_bHat1         (0.0),
  m_muHat1        (1.0),

  m_aAoA1         (0.0),
  m_bAoA1         (0.0),
  m_aAoAHat1      (0.0),
  m_bAoAHat1      (0.0),

  //-------------------------------------------------------------------------//
  // Transient Data:                                                         //
  //-------------------------------------------------------------------------//
  m_mode          (FlightMode::UNDEFINED),
  m_ignTime2      (NAN),
  m_fairingSepTime(NAN),
  m_cutOffTime1   (NAN),
  m_ignTime1      (NAN),
  m_maxQ          (0.0),
  m_sepQ          (0.0),
  m_maxLongG      (0.0),

  //-------------------------------------------------------------------------//
  // Logging Params:                                                         //
  //-------------------------------------------------------------------------//
  m_os            (a_os),        // May be NULL
  m_logLevel      (a_log_level)
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  if (!(0.0 < m_K2            && m_K2 < 1.0            &&
        IsPos(m_fullMass2)    && IsPos(m_emptyMass2)   &&
        IsPos(m_propMass2)    && IsPos(m_unSpendable2) &&
        IsPos(m_spendable2)   &&
        IsPos(m_IspVac2)      && IsPos(m_burnRateI2)   &&
        IsPos(m_thrustVacI2)  && m_thrustMult2 > 0.0   &&
        IsPos(m_T2)           &&
        0.0 <= m_minThrtL2    && m_minThrtL2  <= 1.0   &&
        !IsNeg(m_maxAoA2)     &&
        //
        0.0 < m_K1            && m_K1 < 1.0            &&
        IsPos(m_fullMass1)    && IsPos(m_emptyMass1)   &&
        IsPos(m_propMass1)    && IsPos(m_unSpendable1) &&
        IsPos(m_spendable1)   &&
        IsPos(m_IspVac1)      && IsPos(m_IspSL1)       &&
        m_IspSL1 < m_IspVac1  && IsPos(m_burnRateI1)   &&
        IsPos(m_thrustVacI1)  && m_thrustMult1 > 0.0   &&
        IsPos(m_T1)           &&
        0.0 <= m_minThrtL1    && m_minThrtL1  <= 1.0   &&
        !IsNeg(m_maxAoA1)     &&
        //
        IsPos(m_maxStartMass) && IsPos(m_fairingMass)  &&
        IsPos(m_crosS)))
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
// Returns (RunRC, FinalH, FlightTime, ActStartMass):
//
Ascent2::RunRes Ascent2::Run()
{
  //-------------------------------------------------------------------------//
  // Run the integration BACKWARDS from the orbital insertion point          //
  //-------------------------------------------------------------------------//
  // (@ t=0):
  // Angular LV velocity at the orbital insertion point:
  AngVel omega0 = 1.0_rad * m_Vins / m_Rins;

  // The initial Radial Velocity is 0, and the final delta of Spent Mass
  // is also 0, and the initial Polar Angle is 0:
  StateV s0 = std::make_tuple(m_Rins, VelK(0.0), omega0, 0.0_kg, 0.0_rad);

  // The RHS and the Call-Back Lambdas:
  auto rhs =
    [this](StateV const&  a_s, Time a_t) -> DStateV
    { return this->ODERHS(a_s, a_t); };

  auto cb  =
    [this](StateV* a_s,   Time a_t) -> bool
    { return this->ODECB (a_s, a_t); };

  // FlightMode ctl (modes are switched based on the SpentPropMass):
  m_mode           = FlightMode::Burn2;
  m_ignTime2       = Time(NAN);  // Not known yet
  m_cutOffTime1    = Time(NAN);  // ditto
  m_fairingSepTime = Time(NAN);  // ditto
  m_ignTime1       = Time(NAN);  // ditto

  // The ascent maximum duration (which is certainly enough) is 1 hour:
  constexpr Time t0    = 0.0_sec;
  constexpr Time tMin  = -3600.0_sec;

  // NB: All necessary exception handling is provided inside "RKF5":
  Time tEnd;
  try
  {
    //-----------------------------------------------------------------------//
    // Actually Run the Integrator!                                          //
    //-----------------------------------------------------------------------//
    tEnd =
      RKF5(&s0, t0, tMin, rhs,
           -ODEInitStep, -ODEMaxStep, ODERelPrec, &cb, m_os);
    assert(tEnd >=  tMin);
  }
  catch (NearSingularityExn const& ns)
  {
    //-----------------------------------------------------------------------//
    // We have reached a vicinity of the singular point:                     //
    //-----------------------------------------------------------------------//
    // This is a NORMAL (and moreover, a desirable) outcome. Compute a more
    // precise singular point position:
    //
    auto sing = LocateSingularPoint(ns);

    if (bool(sing))
    {
      StateV const& singS = sing.value().first;
      Time          singT = sing.value().second;
      LenK          singH = std::get<0>(singS) - R;
      assert(!IsNeg(singH));
      Mass          singM = LVMass(singS, singT);

      return RunRes{RunRC::Singularity, singT, singH, VelK(0.0), singM,
                    m_maxQ, m_sepQ, m_maxLongG};
    }
    else
      // Although we got a "NearSingularityExn", we could not determine the
      // precise location of the singular point. Yet the constraints are
      // returned:
      return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN),
                    m_maxQ, m_sepQ, m_maxLongG};
  }
  catch (std::exception const& exn)
  {
    //-----------------------------------------------------------------------//
    // Any other "standard" (or other) exceptions: Ascent unsuccessful:      //
    //-----------------------------------------------------------------------//
    if (m_os != nullptr)
      *m_os  << "# Ascent2::Run: Exception: " << exn.what() << std::endl;

    return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN),
                  m_maxQ, m_sepQ, m_maxLongG};
  }
  catch (...)
  {
    // Any other exception:
    if (m_os != nullptr)
      *m_os  << "# Ascent2::Run: UnKnown Exception" << std::endl;

    return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN),
                  m_maxQ, m_sepQ, m_maxLongG};
  }

  //-------------------------------------------------------------------------//
  // Integration has run to completion, but not to the singular point:       //
  //-------------------------------------------------------------------------//
  // Check the final state and mode:
  //
  LenK   rEnd     = std::get<0>(s0);
  LenK   hEnd     = rEnd - R;
  VelK   VrEnd    = std::get<1>(s0);
  AngVel omegaEnd = std::get<2>(s0);
  VelK   VEnd     = SqRt(Sqr(VrEnd) + Sqr(rEnd * omegaEnd / 1.0_rad));
  Mass   mEnd     = LVMass(s0, tEnd);

  if (m_mode == FlightMode::UNDEFINED)
  {
    // We have run out of propellant while, presumably, still @ hEnd > 0;
    // if we got h <= 0, the neg value should be very small:
    if (!IsPos(hEnd))
    {
      if (m_os != nullptr && m_logLevel >= 1)
        *m_os << "# Ascent2::Run: Got mode=UNDEFINED but h="
              << hEnd.Magnitude() << " km" << std::endl;
      hEnd = 0.0_km;
    }
    return RunRes{RunRC::FlameOut, tEnd, hEnd, VEnd, mEnd,
                  m_maxQ,  m_sepQ, m_maxLongG};
  }
  else
  {
    // Then we should get hEnd==0 (because the only other possibility is
    // that the integration ran until tMin, which is extremely unlikely):
    if (hEnd > 0.5_km)
    {
      if (m_os != nullptr && m_logLevel >= 1)
        *m_os << "# Ascent2::Run: Expected h=0 but got h="
              << hEnd.Magnitude() << " km" << std::endl;
    }
    return RunRes{RunRC::ZeroH, tEnd, 0.0_km, VEnd, mEnd,
                  m_maxQ, m_sepQ, m_maxLongG};
  }
}

//===========================================================================//
// The Mach Number:                                                          //
//===========================================================================//
namespace
{
  inline double Mach(EAM::AtmConds const& a_atm, VelK a_v)
  {
    assert(!IsNeg(a_v));

    // The Speed of Sound:
    Vel A = std::get<3>(a_atm);
    return
      IsZero(A)
      ? Inf<double>                 // We are above the atmosphere
      : double(To_Len_m(a_v) / A);  // Genertic case
  }
}

//===========================================================================//
// "NonGravForces":                                                          //
//===========================================================================//
// Common Helper for "ODERHS" and "ODECB":
//
void Ascent2::NonGravForces
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
  Angle aoa    = AoA(a_t, psi);
  assert(!IsNeg(aoa));
  double cosA  = Cos(aoa);
  double sinA  = Sin(aoa);

  // AeroDynamic Drag and Lift:
  auto [atm, drag, lift]  = AeroDynForces(a_r, V, aoa);
  assert(!(IsNeg(drag) || IsNeg(lift)));

  // Propellant Burn Rate (>= 0) and Thrust:
  MassRate burnRate = PropBurnRate(a_t);
  assert(!IsNeg(burnRate));

  Pressure p        = std::get<0> (atm);
  ForceK   thrust   = Thrust(burnRate, p);
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
Ascent2::DStateV Ascent2::ODERHS(StateV const& a_s, Time a_t)
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
  Mass   m  = LVMass(a_s, a_t);
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
// ODE CallBack (invoked after the completion of each RKF5 step):            //
//===========================================================================//
// Here FlightMode switching occurs, so this method is non-"const":
//
bool Ascent2::ODECB(StateV* a_s, Time a_t)
{
  assert(a_s != nullptr && !IsPos(a_t));
  LenK   r             = std::get<0>(*a_s);
  LenK   h             = r - R;             // Altitude
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  Mass   spentPropMass = std::get<3>(*a_s);
  Angle  phi           = std::get<4>(*a_s);

  // NB: "spentPropMass" decreases over time, so increases in Bwd time.
  // It is 0 at Orbit Insertion Time (t=0):
  assert(IsPos(r) && !IsNeg(spentPropMass));

  //-------------------------------------------------------------------------//
  // Similar to "ODERHS", check if we are approaching the singularity:       //
  //-------------------------------------------------------------------------//
  VelK   Vhor          = r * omega / 1.0_rad;

  if (Vhor < SingVhor)
  {
    omega              = AngVel(0.0);
    std::get<2>(*a_s)  = AngVel(0.0);
    Vhor               = VelK(0.0);
  }
  if (IsZero(omega) && Vr < SingVr)
    throw NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

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
  double   M           = Mach(atm, V);
  Density  rho         = std::get<1>(atm);
  auto     V2          = To_Len_m(Sqr(V));
  Pressure Q           = 0.5 * rho * V2;

  // Output at the beginning:
  if (IsZero(a_t) && m_os != nullptr && m_logLevel >= 1)
  {
    assert(m_mode == FlightMode::Burn2);
    *m_os << "# t=0 sec, h=" << h.Magnitude()
          << " km, V="       << V.Magnitude() << " km/sec, Mass="
          << m.Magnitude()   << " kg"         << std::endl;
  }
  //-------------------------------------------------------------------------//
  // Switching Burn2 -> Gap:                                                 //
  //-------------------------------------------------------------------------//
  // Occurs according to "spentPropMass": when all Stage2 propellants (except
  // the Remnants) is spent:
  //
  if (m_mode == FlightMode::Burn2 && spentPropMass >= m_spendable2)
  {
    assert(!IsFinite(m_ignTime2));
    m_ignTime2    = a_t;
    assert(IsNeg(m_ignTime2));
    m_mode        = FlightMode::Gap;
    // Furthermore, since the duration of the Ballistic Gap is known, at
    // this point we already know the Stage1 Cut-Off time:
    // NB:  HERE "m_TGap" is used:
    assert(!IsNeg(m_TGap) && !IsFinite(m_cutOffTime1));
    m_cutOffTime1 = a_t - m_TGap;
    assert(m_cutOffTime1 <= m_ignTime2);

    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="              << psi_deg.Magnitude()
            << " deg, m=" << m.Magnitude()  << " kg: "
               "Stage2 Ignition, Ballistic Gap Ends"      << std::endl;
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
      m_mode   = FlightMode::Burn1;

      // XXX: We assume that Stage1 separation occures at the cut-off moment,
      // so record the corresp Mach number:
      m_sepQ   = Q;

      if (m_os != nullptr && m_logLevel >= 1)
      {
        // Re-calculate the Mass in Burn1 for the output (now with Stage1):
        m      = LVMass(*a_s, a_t);
        *m_os << "# t="     << a_t.Magnitude() << " sec, h=" << h.Magnitude()
              << " km, L="  << L.Magnitude()   << " km, V="  << V.Magnitude()
              << " km/sec, psi="               << psi_deg.Magnitude()
              << " deg, m=" << m.Magnitude()   << " kg, M= " << M
              << ", Q="     << Q.Magnitude()
              << ": Stage1 Cut-Off and Separation,  Ballistic Gap Starts"
              << std::endl;
      }
    }
  }
  //-------------------------------------------------------------------------//
  // Switching Burn1 -> UNDEFINED:                                           //
  //-------------------------------------------------------------------------//
  // When all Stage1 propellants (except the Remnants) is spent:
  //
  if (m_mode == FlightMode::Burn1  &&
      spentPropMass >= m_spendable2 + m_spendable1)
  {
    assert(!IsFinite(m_ignTime1));
    m_ignTime1 = a_t;
    m_mode     = FlightMode::UNDEFINED;
    assert(m_ignTime1 < m_cutOffTime1 && m_cutOffTime1 <= m_ignTime2 &&
           IsNeg(m_ignTime2));

    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# t="     << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L="  << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="               << psi_deg.Magnitude()
            << " deg, m=" << m.Magnitude()   << " kg: "
               "Stage1 Ignition"             << std::endl;
  }

  //-------------------------------------------------------------------------//
  // Monitor the Constraints:                                                //
  //-------------------------------------------------------------------------//
  // Dynamic Pressure:
  m_maxQ            = std::max(Q, m_maxQ);

  // LongG:
  // This is a projection of "ngAcc" to the main LV axis (XXX: verify whether
  // this expr is also correct in case of Falling Back to Earth).
  // In the (r, normal_to_r) frame, the longitudinal LV axis is given by the
  // (psi + aoa) angle:
  AccK   longAcc = lvAxis[0] * ngAcc[0] + lvAxis[1] * ngAcc[1];
  double longG   = double  (longAcc / g0K);
  m_maxLongG     = std::max(longG,  m_maxLongG);

  //-------------------------------------------------------------------------//
  // In any mode, detect the Fairing Separation Condition:                   //
  //-------------------------------------------------------------------------//
  // In the reverse time, it's when the Dynamic Pressure becomes HIGHER
  // that the threshold:
  if (!IsFinite(m_fairingSepTime) && Q >= FairingSepCond)
  {
    m_fairingSepTime = a_t;

    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# t="     << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L="  << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="               << psi_deg.Magnitude()
            << " deg, m=" << m.Magnitude()   << " kg, M="  << M
            << ", Q="     << Q.Magnitude()   << ": Fairing Separation"
            << std::endl;
  }
  //---------------------------------------------------------------------//
  // Stopping Conds:                                                     //
  //---------------------------------------------------------------------//
  bool cont = true;

  // Do not continue if:
  // (*) we are on the surface (stopping by "h");
  // (*) we are in the UNDEFINED mode (stopping by Ignition Time);
  if (!IsPos(h))
  {
    cont = false;

    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# STOP: t="   << a_t.Magnitude()  << " sec: H=0, Vr="
            << Vr.Magnitude() << " km/sec, Vhor=" << Vhor.Magnitude()
            << " km/sec, m="  << m.Magnitude()    << " kg, L="
            << L.Magnitude()  << " km, psi="      << psi_deg.Magnitude()
            << " deg, "       << ToString(m_mode) << std::endl;
  }
  if (m_mode == FlightMode::UNDEFINED)
  {
    cont = false;

    if (m_os != nullptr && m_logLevel >= 1)
      *m_os << "# STOP: t="   << a_t.Magnitude()  << " sec: UNDEF, Vr="
            << Vr.Magnitude() << " km/sec, Vhor=" << Vhor.Magnitude()
            << " km/sec, m="  << m.Magnitude()    << " kg, L="
            << L.Magnitude()  << " km, psi="      << psi_deg.Magnitude()
            << " deg, "       << ToString(m_mode) << std::endl;
  }

  // XXX: OUTPUT (with a 100 msec step, or if we are going to stop now):
  if (m_os  != nullptr && m_logLevel >= 3 &&
     (!cont || int(Round(double(a_t / 0.001_sec))) % 100 == 0))
  {
    // Thrust is more conveniently reported in kgf:
    auto tkg = thrust / g0K;

    *m_os << a_t.Magnitude()     << '\t' << h.Magnitude()        << '\t'
          << L.Magnitude()       << '\t' << Vr.Magnitude()       << '\t'
          << Vhor.Magnitude()    << '\t' << V.Magnitude()        << '\t'
          << psi_deg.Magnitude() << '\t' << aoa_deg.Magnitude()  << '\t'
          << m.Magnitude()       << '\t' << ToString(m_mode)     << '\t'
          << tkg.Magnitude()     << '\t' << burnRate.Magnitude() << '\t'
          << Q.Magnitude()       << '\t' << M                    << '\t'
          << longG               << std::endl;
  }
  return cont;
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
Ascent2::AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const
{
  auto atm  = EAM::GetAtmConds(std::max(a_r - R, 0.0_km));
  double M  = Mach(atm, a_v);
  if (!IsFinite(M))
    return std::make_tuple(atm, ForceK(0.0), ForceK(0.0));

  // Generic Case:
  assert(M >= 0.0);

  // The Aerodynamic Force Main Term:
  ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * m_crosS;

  // The Drag and Lift Coeffs, using the default model:
  double cD = LVAeroDyn::cD(M, a_AoA);
  double cL = LVAeroDyn::cL(M, a_AoA);

  return std::make_tuple(atm, cD * F, cL * F);
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
std::optional<std::pair<Ascent2::StateV, Time>>
Ascent2::LocateSingularPoint(NearSingularityExn const& a_nse)
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
  MassRate burnRate1 = PropBurnRate(t1);
  assert(!IsNeg(burnRate1));

  ForceK thrust1 = Thrust(burnRate1, std::get<0>(atm));
  assert(!IsNeg(thrust1) && IsZero(thrust1) == IsZero(burnRate1));

  // The Curr LV Mass (via a synthetic "s1"):
  StateV s1    { r1, Vr1, AngVel(0.0), a_nse.m_spentPropMass, a_nse.m_phi };
  Mass   m1  = LVMass(s1, t1);
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
    *m_os << "# SingularPoint Located: t1="     << t1.Magnitude()
          << " sec, tau="   << tau.Magnitude()  << " sec, tS="
          << tS.Magnitude() << " sec, m1="      << m1.Magnitude()
          << " kg, mS="     << (m1 + burnRate1 * tau).Magnitude()
          << " kg, mSalt="  << LVMass(singS, tS).Magnitude()
          << " kg, hS="     << (rS - R).Magnitude()
          << " km, LS="     << LS      .Magnitude()
          << " km, "        << ToString(m_mode) << std::endl;

  return std::make_optional(std::make_pair(singS, tS));
}

//===========================================================================//
// Propellant Burn Rate: May be variable:                                    //
//===========================================================================//
MassRate Ascent2::PropBurnRate(Time a_t) const
{
  switch (m_mode)
  {
  case FlightMode::Burn1:
  {
    // Arg: "tau" is the time since Stage1 ignition  (so 0 <= t <= T1); but the
    // latter is not yet reached, so calculate it using "T1";
    Time     tIgn1 = m_cutOffTime1  - m_T1;
    Time     tau   = std::min(std::max(a_t  - tIgn1,    0.0_sec), m_T1);
    MassRate res   = std::max(m_burnRateI1  + (m_bMu1 + m_aMu1 * tau) * tau,
               MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));
    return res;
  }
  case FlightMode::Gap:
    return MassRate(0.0);

  case FlightMode::Burn2:
  {
    // Arg: "tau" is the time since Stage2 ignition  (so 0 <= t <= T2); but the
    // latter is not yet reached, so calculate it using "T2":
    Time     tIgn2 = - m_T2;
    Time     tau   = std::min(std::max(a_t  - tIgn2,    0.0_sec), m_T2);
    MassRate res   = std::max(m_burnRateI2  + (m_bMu2 + m_aMu2 * tau) * tau,
                              MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));
    return res;
  }
  default:
    return MassRate(0.0);
  }
}

//===========================================================================//
// Angle-of-Attack:                                                          //
//===========================================================================//
Angle Ascent2::AoA(Time a_t, Angle a_psi) const
{
  Angle  aoa = 0.0_rad;

  switch (m_mode)
  {
  case FlightMode::Burn1:
  {
    // Arg: "tau" which is the time since Stage1 ignition (where the latter is
    // calculated via "m_T1", because the actual event-based "m_ignTime1"   is
    // not reached yet), so AoA1(tau=0)=0, ie @ launch;   0 <= tau <= T1:
    Time tIgn1 = m_cutOffTime1  - m_T1;
    Time tau   = std::min(std::max(a_t - tIgn1,    0.0_sec), m_T1);
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
    return 0.0_rad;

  case FlightMode::Burn2:
  {
    // Arg: "t" is time from Stage2 cut-off (same as orbital insertion time),
    // so it is just the main time "a_t": -T2 <= t <= 0; enforce those constr-
    // aints. NB: Formulas similar to those for Stage2, but with the INVERTED
    // sign of "b"; AoA2(t=0)=0, ie @ orbital insetrion:
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

  default:
    return 0.0_rad;
  }
  // Check that (AoA + TrajIncl <= Pi/2), ie, the AoA does not point the thrust
  // backwards  (XXX: there is no similar constraint if we are descending):
  if (aoa + a_psi > PI_2)
    aoa = PI_2 - a_psi;
  return aoa;
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
      Time     Isp = x   * m_IspSL1 + (1.0 - x) * m_IspVac1;
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
// Current LV Mass (incl the PayLoad):                                       //
//===========================================================================//
Mass Ascent2::LVMass(StateV const& a_s, Time a_t) const
{
  // The mass of Propellants spent between "a_t" < 0 and the Orbital Inser-
  // tion instant (t=0):
  Mass spentPropMass = std::get<3>(a_s);
  assert(!IsNeg(spentPropMass));

  // Initialise "m" to the mass @ Orbital Insertion (incl the unspent Stage2
  // Propellant):
  Mass m = m_emptyMass2 + m_unSpendable2 + m_payLoadMass;

  // In any case, ADD the SpentPropMass (between "a_t" and Orbital Insertion):
  m += spentPropMass;

  // If Stage1 has NOT separated yet, add its Empty and PropRemnants:
  if (m_mode != FlightMode::Burn2 && m_mode != FlightMode::Gap)
    m += m_emptyMass1 + m_unSpendable1;

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
  if (m_os == nullptr)
    return;

  // Here "t" is the time from Stage2 cut-off (so t=-T2..0):
  *m_os << "# T2   := " << m_T2.Magnitude() << ';'    << std::endl;
  *m_os << "# 0  <= tau <= T2"                        << std::endl;
  *m_os << "#-T2 <= t   <= 0"                         << std::endl;
  *m_os << "# AoA2 := t * ("   << m_aAoA2.Magnitude() << " * t - ("
        << m_bAoA2.Magnitude() << ")); "              << std::endl;
  *m_os << "# mu2  := " << m_burnRateI2.Magnitude()   << " + ("
        << m_bMu2.Magnitude()  << ") * tau + ("
        << m_aMu2.Magnitude()  << ") * tau^2;"        << std::endl;

  // Here "t" is the time from Stage1 cut-off (so t=-T1..0), and "tau" is
  // the time since Stage1 ignition (NB: the (-b) coeff at "tau" for AoA),
  // so tau=0..T1:
  *m_os << "# T1   := " << m_T1.Magnitude() << ';'    << std::endl;
  *m_os << "# 0  <= tau <= T1"                        << std::endl;
  *m_os << "# AoA1 := tau * (" << m_aAoA1.Magnitude() << " * tau + ("
        << m_bAoA1.Magnitude() << ")); "              << std::endl;
  *m_os << "# mu1  := " << m_burnRateI1.Magnitude()   << " + ("
        << m_bMu1.Magnitude()  << ") * tau + ("
        << m_aMu1.Magnitude()  << ") * tau^2;"        << std::endl;
}
}
// End namespace SpaceBallistics
