// vim:ts=2:et
//===========================================================================//
//                          "Src/Missions/Ascent2.cpp":                      //
//                   Ascent-to-Orbit for a "Model" 2-Stage LV                //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"
#include "SpaceBallistics/PhysEffects/LVAeroDyn.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"
#include <nlopt.hpp>

namespace SpaceBallistics
{
//===========================================================================//
// Static Cache:                                                             //
//===========================================================================//
std::unordered_map<Ascent2::AscCtls, Ascent2::RunRes> Ascent2::s_Cache;

//===========================================================================//
// "Ascent2": Non-Default Ctor:                                              //
//===========================================================================//
Ascent2::Ascent2
(
  // LV Params:
  double          a_alpha1,      // FullMass1 / FullMass2
  ForceK          a_thrust2_vac,
  ForceK          a_thrust1_vac,

  // Mission Params:
  Mass            a_payload_mass,
  LenK            a_h_perigee,
  LenK            a_h_apogee,
  Angle_deg       a_incl,
  Angle_deg       a_launch_lat,
  AscCtls const&  a_ctls,        // Flight Ctls

  // Logging Params:
  std::ostream*   a_os,
  int             a_log_level
)
: //-------------------------------------------------------------------------//
  // LV Params:                                                              //
  //-------------------------------------------------------------------------//
  // Stage2 Params:
  m_fullMass2     ((MaxStartMass - FairingMass - a_payload_mass) /
                   (1.0 + a_alpha1)),
  m_emptyMass2    (m_fullMass2   * (1.0 - K2)),
  m_propMass2     (m_fullMass2   * K2),
  m_thrustVac2    (a_thrust2_vac),
  m_burnRateI2    (m_thrustVac2  / (IspVac2  * g0K)),

  // Stage1 Params:
  m_fullMass1     ((MaxStartMass - FairingMass - a_payload_mass) /
                   (1.0 + a_alpha1) * a_alpha1),
  m_emptyMass1    (m_fullMass1  * (1.0 - K1)),
  m_propMass1     (m_fullMass1  * K1),
  m_thrustVac1    (a_thrust1_vac),
  m_burnRateI1    (m_thrustVac1  / (IspVac1  * g0K)),

  //-------------------------------------------------------------------------//
  // Mission Params:                                                         //
  //-------------------------------------------------------------------------//
  m_payLoadMass   (a_payload_mass),
  m_Rins          (),            // Not yet...
  m_Vins          (),            //

  //-------------------------------------------------------------------------//
  // Flight Ctl Params:                                                      //
  //-------------------------------------------------------------------------//
  m_TGap          (),            // 0 by default

  // BurnRate Ctls and AoA Ctls. NB: "m_T{1,2}" are initialised to Actual
  // BurnTimes (taking the Remnants into account):
  m_T2            (m_propMass2 * (1.0 - PropRem2) / m_burnRateI2),
  m_aMu2          (MassT3(0.0)),
  m_bMu2          (MassT2(0.0)),

  m_T1            (m_propMass1 * (1.0 - PropRem1) / m_burnRateI1),
  m_aMu1          (MassT3(0.0)),
  m_bMu1          (MassT2(0.0)),

  m_aAoA2         (AngAcc(0.0)),
  m_bAoA2         (AngVel(0.0)),
  m_aAoA1         (AngAcc(0.0)),
  m_bAoA1         (AngVel(0.0)),

  //-------------------------------------------------------------------------//
  // Transient Data:                                                         //
  //-------------------------------------------------------------------------//
  m_mode          (FlightMode::UNDEFINED),
  m_ignTime2      (Time(NAN)),
  m_fairingSepTime(Time(NAN)),
  m_cutOffTime1   (Time(NAN)),
  m_ignTime1      (Time(NAN)),

  //------------------------------------------------------------------------//
  // Logging Params:                                                        //
  //------------------------------------------------------------------------//
  m_os            (a_os),        // May be NULL
  m_logLevel      (a_log_level)
{
  //-------------------------------------------------------------------------//
  // Check the Params:                                                       //
  //-------------------------------------------------------------------------//
  if ( IsNeg(m_payLoadMass) || !IsPos(m_fullMass2) || !IsPos(m_thrustVac2)  ||
      !IsPos(m_fullMass1)   || !IsPos(m_thrustVac1))
    throw std::invalid_argument("Ascent2::Ctor: Invalid LV param(s)");

  // Therefore:
  assert(IsPos(m_emptyMass2) && IsPos(m_propMass2) && IsPos(m_burnRateI2) &&
         IsPos(m_emptyMass1) && IsPos(m_propMass1) && IsPos(m_burnRateI1));
  if (a_h_perigee  <  100.0_km || a_h_perigee  > a_h_apogee ||
      a_incl       <  0.0_deg  || a_incl       > 180.0_deg  ||
      // NB: We assume launch from the Northern Hemisphere, and not from the
      // North Pole:
      a_launch_lat <   0.0_deg || a_launch_lat >= 90.0_deg)
    throw std::invalid_argument
          ("Ascent2::Ctor: Invalid Orbit/Launch  Param(s)");

  //-------------------------------------------------------------------------//
  // Orbital and Launch Params:                                              //
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

  //-------------------------------------------------------------------------//
  // Finally, set the Flight Ctl Params:                                     //
  //-------------------------------------------------------------------------//
  SetCtlParams(a_ctls);
}

//===========================================================================//
// "SetCtlParams":                                                           //
//===========================================================================//
// With "AscCtls" struct:
//
void Ascent2::SetCtlParams(AscCtls const& a_ctls)
{
  SetCtlParams(a_ctls.m_aHat2, a_ctls.m_muHat2,
               a_ctls.m_aAoA2, a_ctls.m_bAoA2,
               a_ctls.m_TGap,
               a_ctls.m_aHat1, a_ctls.m_muHat1,
               a_ctls.m_aAoA1, a_ctls.m_bAoA1);
}

// With individual params: The actual implementation:
//
void Ascent2::SetCtlParams
(
  double a_aHat2, double a_muHat2,  // Stage2 BurnRate
  double a_aAoA2, double a_bAoA2,   // Stage2 AoA
  double a_TGap,                    // Ballistic Gap in sec
  double a_aHat1, double a_muHat1,  // Stage1 BurnRate
  double a_aAoA1, double a_bAoA1    // Stage1 AoA
)
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  if (a_TGap   <  0.0                      ||
      std::fabs(a_aHat2)            > 1.0  ||
      a_muHat2 <= 0.0   || a_muHat2 > 1.0  ||
      a_aAoA2  <  0.0   || a_aAoA2  > 1.0  ||
      a_bAoA2  <  0.0   || a_bAoA2  > 1.0  ||
      std::fabs(a_aHat1)            > 1.0  ||
      a_muHat1 <= 0.0   || a_muHat2 > 1.0  ||
      a_aAoA1  <  0.0   || a_aAoA1  > 1.0  ||
      a_bAoA1  <  0.0   || a_bAoA1  > 1.0)
      throw std::invalid_argument("Ascent2::SetCtlParams: Invalid Param(s)");

  //-------------------------------------------------------------------------//
  // Ballistic Gap:                                                          //
  //-------------------------------------------------------------------------//
  m_TGap = Time(a_TGap);

  //-------------------------------------------------------------------------//
  // BurnRate Coeffs:                                                        //
  //-------------------------------------------------------------------------//
  // Normalise "aHat"s to [-1/3 .. +1/3], then compute the actual limits for
  // "muHat", and the actual "muHat" within those limits:
  // Stage2:
  double aHat2    = a_aHat2 / 3.0;
  double muHatLo2 = (1.0 - aHat2)  / 2.0;
  double muHatUp2 = (aHat2 < 0.0) ? 1.0 + aHat2 : 1.0 - 2.0 * aHat2;
  assert(muHatLo2 <= muHatUp2);
  double muHat2   = a_muHat2 * (muHatUp2 - muHatLo2);
  // XXX: In a very degenerate case, we may get muHat2==0:
  assert(0.0 < muHat2 && muHat2 <= 1.0);
  m_T2            = m_propMass2 / (m_burnRateI2    * muHat2);
  m_aMu2          = 3.0 * m_burnRateI2 / Sqr(m_T2) * aHat2;
  m_bMu2          = 2.0 * (m_propMass2 / Sqr(m_T2) - m_burnRateI2 / m_T2 -
                           m_aMu2 * m_T2 / 3.0);

  // Stage1:
  double aHat1    = a_aHat1 / 3.0;
  double muHatLo1 = (1.0 - aHat1)  / 2.0;
  double muHatUp1 = (aHat1 < 0.0) ? 1.0 + aHat1 : 1.0 - 2.0 * aHat1;
  assert(muHatLo1 <= muHatUp1);
  double muHat1   = a_muHat1 * (muHatUp1 - muHatLo1);
  // XXX: In a very degenerate case, we may get muHat1==0:
  assert(0.0 < muHat1 && muHat1 <= 1.0);
  m_T1            = m_propMass1 / (m_burnRateI1    * muHat1);
  m_aMu1          = 3.0 * m_burnRateI1 / Sqr(m_T1) * aHat1;
  m_bMu1          = 2.0 * (m_propMass1 / Sqr(m_T1) - m_burnRateI1 / m_T1 -
                           m_aMu1 * m_T1 / 3.0);

  //-------------------------------------------------------------------------//
  // AoA Coeffs:                                                             //
  //-------------------------------------------------------------------------//
  // Stage2:
  m_bAoA2        = - 4.0 * MaxAoA2 / m_T2 * a_bAoA2;
  AngAcc aAoALo2 = m_bAoA2 / m_T2;
  AngAcc aAoAUp2 =
    (a_bAoA2 < 0.5)
    ? (MaxAoA2 / m_T2 + m_bAoA2) / m_T2
    : - Sqr(m_bAoA2) / (4.0 * MaxAoA2);
  m_aAoA2 = aAoALo2  * (1.0 - a_aAoA2) + a_aAoA2 * aAoAUp2;

  // Stage1:
  m_bAoA1        = - 4.0 * MaxAoA1 / m_T1 * a_bAoA1;
  AngAcc aAoALo1 = m_bAoA1 / m_T1;
  AngAcc aAoAUp1 =
    (a_bAoA1 < 0.5)
    ? (MaxAoA1 / m_T1 + m_bAoA1) / m_T1
    : - Sqr(m_bAoA1) / (4.0 * MaxAoA1);
  m_aAoA1 = aAoALo1  * (1.0 - a_aAoA1) + a_aAoA1 * aAoAUp1;
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

      return RunRes{RunRC::Singularity, singT, singH, VelK(0.0), singM};
    }
    else
      // Although we got a "NearSingularityExn", we could not determine the
      // precise location of the singular point:
      return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN)};
  }
  catch (std::exception const& exn)
  {
    //-----------------------------------------------------------------------//
    // Any other "standard" (or other) exceptions: Ascent unsuccessful:      //
    //-----------------------------------------------------------------------//
    if (m_os != nullptr)
      *m_os  << "# Ascent2::Run: Exception: " << exn.what() << std::endl;

    return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN)};
  }
  catch (...)
  {
    // Any other exception:
    if (m_os != nullptr)
      *m_os  << "# Ascent2::Run: UnKnown Exception" << std::endl;

    return RunRes{RunRC::Error, Time(NAN), LenK(NAN), VelK(NAN), Mass(NAN)};
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
    return RunRes{RunRC::FlameOut, tEnd, hEnd, VEnd, mEnd};
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
    return RunRes{RunRC::ZeroH, tEnd, 0.0_km, VEnd, mEnd};
  }
}

//===========================================================================//
// ODE RHS for a simple Gravity Turn:                                        //
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
  // Abs Velocity: Since we are NOT near the singularity, it should be
  // sufficiently far away from 0:
  VelK   V  = SqRt(Sqr(Vr) + Sqr(Vhor));
  assert(IsPos(V));

  // The Angle-of-Attack:
  Angle  aoa  = AoA(a_t);
  assert(!IsNeg(aoa));
  double cosA = Cos(aoa);
  double sinA = Sin(aoa);

  // AeroDynamic Drag and Lift:
  auto [atm, drag, lift]  = AeroDynForces(r, V, aoa);
  assert(!(IsNeg(drag) || IsNeg(lift)));

  // Propellant Burn Rate (>= 0) and Thrust:
  MassRate burnRate = PropBurnRate(a_t);
  assert(!IsNeg(burnRate));

  ForceK   thrust   = Thrust(burnRate, std::get<0>(atm));
  assert(!IsNeg(thrust) && IsZero(thrust) == IsZero(burnRate));

  // The curr mass:
  Mass   m  = LVMass(a_s, a_t);
  assert(IsPos(m));

  // "u" unit vector: normal to the velocity (in the "up" direction)
  // "v" unit vector: in the velocity direction
  // Components of "u" in the (radius-vector, normal-to-radius-vector)
  // frame:
  double v[2] { double(Vr / V), double(r * omega / V) };
  double u[2] { v[1],           -v[0]                 };

  // Thrust components (in the same frame):
  ForceK T[2] { (v[0] * cosA + u[0] * sinA) * thrust,
                (v[1] * cosA + u[1] * sinA) * thrust };

  // AeroDynamic Drag (-v) and Lift (+u) components:
  ForceK D[2] { -v[0] * drag, -v[1] * drag };
  ForceK L[2] {  u[0] * lift,  u[1] * lift };

  // Acceleration components due to the above forces:
  AccK acc[2] { (T[0] + D[0] + L[0]) / m,
                (T[1] + D[1] + L[1]) / m };

  // The RHS components:
  AccK r2Dot =
      r * Sqr(omega / 1.0_rad)   // "Kinematic"     term
    - K / Sqr(r)                 // Gravitational   term
    + acc[0];                    // Thrust-, Drag-, Lift-indiced accs

  AngAcc omegaDot =
    (
      - 2.0 * Vr * omega         // "Kinematic"     term
      + acc[1]   * 1.0_rad       // Thrust-, Drag-, Lift-indiced accs
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
  LenK   h             = r - R;
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  Mass   spentPropMass = std::get<3>(*a_s);
  Mass   m             = LVMass(*a_s, a_t); // Using the "old" Mode yet!
  Angle  phi           = std::get<4>(*a_s);

  //-------------------------------------------------------------------------//
  // Similar to "ODERHS", check if we are approaching the singularity:       //
  //-------------------------------------------------------------------------//
  VelK   Vhor          = r * omega / 1.0_rad;
  auto   V2K           = Sqr(Vr) + Sqr(Vhor);
  VelK   V             = SqRt(V2K);
  if (Vhor < SingVhor)
  {
    omega             = AngVel(0.0);
    std::get<2>(*a_s) = AngVel(0.0);
    Vhor              = VelK(0.0);
  }
  if (IsZero(omega) && Vr < SingVr)
    throw NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  //-------------------------------------------------------------------------//
  // Generic Case:                                                           //
  //-------------------------------------------------------------------------//
  // NB: "spentPropMass" decreases over time, so increases in Bwd time.
  // It is 0 at Orbit Insertion Time (t=0):
  assert(IsPos(r) && !IsNeg(spentPropMass) && IsPos(m) && IsPos(V));

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
  if (m_mode == FlightMode::Burn2  &&
      spentPropMass >= m_propMass2 * (1.0 - PropRem2))
  {
    m_ignTime2    = a_t;
    assert(IsNeg(m_ignTime2));
    m_mode        = FlightMode::Gap;
    // Furthermore, since the duration of the Ballistic Gap is known, at
    // this point we already know the Stage1 Cut-Off time:
    // NB: HERE "m_TGap" is used:
    assert(!IsNeg(m_TGap));
    m_cutOffTime1 = a_t - m_TGap;
    assert(m_cutOffTime1 <= m_ignTime2);

    if (m_os != nullptr && m_logLevel >= 1)
    {
      // It might useful to output the Pitch (= Trajectory Inclination in
      // this case) Angle, and DownRange distance (along the Earth surface):
      double psi = ATan2(Vr, Vhor);
      LenK   L   = R * double(phi);    // Will be < 0

      *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="              << psi        << ", m="
            << m.Magnitude()                << " kg: "
               "Stage2 Ignition, Ballistic Gap Ends"      << std::endl;
    }
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

      if (m_os != nullptr && m_logLevel >= 1)
      {
        // Re-calculate the Mass in Burn1 for the output (with Stage1);
        // again, output "L" and "psi":
        m          = LVMass(*a_s, a_t);
        double psi = ATan2(Vr, Vhor);
        LenK   L   = R * double(phi);  // Will be < 0

        *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
              << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
              << " km/sec, psi="              << psi        << ", m="
              << m.Magnitude()                << " kg: "
                   "Stage1 Cut-Off,  Ballistic Gap Starts"  << std::endl;
      }
    }
  }
  //-------------------------------------------------------------------------//
  // Switching Burn1 -> UNDEFINED:                                           //
  //-------------------------------------------------------------------------//
  // When all Stage1 propellants (except the Remnants) is spent:
  //
  if (m_mode == FlightMode::Burn1 &&
      spentPropMass >= m_propMass2 * (1.0 - PropRem2) +
                       m_propMass1 * (1.0 - PropRem1))
  {
    m_ignTime1 = a_t;
    m_mode     = FlightMode::UNDEFINED;
    assert(m_ignTime1 < m_cutOffTime1 && m_cutOffTime1 <= m_ignTime2 &&
           IsNeg(m_ignTime2));

    if (m_os != nullptr && m_logLevel >= 1)
    {
      // Again, output "L" and "psi":
      double psi = ATan2(Vr, Vhor);
      LenK   L   = R * double(phi);  // Will be < 0

      *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="              << psi        << ", m="
            << m.Magnitude()                << " kg: "
               "Stage1 Ignition"            << std::endl;
    }
  }
  //-------------------------------------------------------------------------//
  // In any mode, detect the Fairing Separation Condition:                   //
  //-------------------------------------------------------------------------//
  // In the reverse time, it's when the Dynamic Pressure becomes HIGHER
  // that the threshold:
  auto     atm = EAM::GetAtmConds(std::max(h, 0.0_km));
  Density  rho = std::get<1>(atm);
  auto     V2  = To_Len_m(V2K);
  Pressure Q   = 0.5 * rho * V2;

  if (IsZero(m_fairingSepTime) && Q >= FairingSepCond)
  {
    m_fairingSepTime = a_t;

    if (m_os != nullptr && m_logLevel >= 1)
    {
      // Again, output "L" and "psi":
      double psi = ATan2(Vr, Vhor);
      LenK   L   = R * double(phi);  // Will be < 0

      *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
            << " km, L=" << L.Magnitude()   << " km, V="  << V.Magnitude()
            << " km/sec, psi="              << psi        << ", m="
            << m.Magnitude()                << " kg: "
               "Fairing Separation"         << std::endl;
    }
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
    {
      // Again, output "L" and "psi":
      double psi = ATan2(Vr, Vhor);
      LenK   L   = R * double(phi);  // Will be < 0

      *m_os << "# t=" << a_t.Magnitude()          << " sec: H=0, Vr="
            << Vr.Magnitude() << " km/sec, Vhor=" << Vhor.Magnitude()
            << " km/sec, m="  << m.Magnitude()    << " kg, L="
            << L.Magnitude()  << " km, psi="      << psi
            << ", "           << ToString(m_mode) << std::endl;
    }
  }
  if (m_mode == FlightMode::UNDEFINED)
  {
    cont = false;

    if (m_os != nullptr && m_logLevel >= 1)
    {
      // Again, output "L" and "psi":
      double psi = ATan2(Vr, Vhor);
      LenK   L   = R * double(phi);  // Will be < 0

      *m_os << "# t=" << a_t.Magnitude()            << " sec: UNDEF, Vr="
            << Vr.Magnitude()   << " km/sec, Vhor=" << Vhor.Magnitude()
            << " km/sec, m="    << m.Magnitude()    << " kg, L="
            << L.Magnitude()    << " km, psi="      << psi
            << ", "             << ToString(m_mode) << std::endl;
    }
  }

  // XXX: Output (with a 100 msec step, or if we are going to stop now):
  if (m_os  != nullptr && m_logLevel >= 3 &&
     (!cont || int(Round(double(a_t / 0.001_sec))) % 100 == 0))
  {
    // Thrust and BurnRate are for info only:
    MassRate burnRate = PropBurnRate(a_t);
    Pressure p        = std::get<0>(atm);
    auto     thrust   = Thrust(burnRate, p) / g0K;

    *m_os << a_t.Magnitude()      << '\t' << h.Magnitude()      << '\t'
          << Vr.Magnitude()       << '\t' << Vhor.Magnitude()   << '\t'
          << V.Magnitude()        << '\t' << m.Magnitude()      << '\t'
          << ToString(m_mode)     << '\t' << thrust.Magnitude() << '\t'
          << burnRate.Magnitude() << '\t' << Q.Magnitude()
          << std::endl;
  }
  return cont;
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
Ascent2::AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA)
{
  // Curr altitude (possible slightly negative vals are rounded to 0):
  auto atm  = EAM::GetAtmConds(a_r - R);

  // The Speed of Sound:
  Vel  A    = std::get<3>(atm);
  if (IsZero(A))
    // We are above the atmosphere:
    return std::make_tuple(atm, ForceK(0.0), ForceK(0.0));

  // The Mach Number:
  double M  = double(To_Len_m(a_v) / A);
  assert(M >= 0.0);

  // The Aerodynamic Force Main Term:
  ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * CroS;

  // The Drag and Lift Coeffs:
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
  // Mass is CONSTANT at this final short interval:
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
    // Then the LV is falling back to the pad, and the singular point is not
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
    // Arg: "t" is the time from Stage1 cut-off (so -T1 <= t <= 0):
    Time     t   = std::max(std::min(a_t - m_cutOffTime1, 0.0_sec), -m_T1);
    MassRate res = std::max(m_burnRateI1 + (m_bMu1 + m_aMu1 * t) * t, 
                            MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));
    return res;
  }
  case FlightMode::Gap:
    return MassRate(0.0);

  case FlightMode::Burn2:
  {
    // Arg: "t" is time from Stage2 cut-off (same as orbital insertion
    // time), so it is just the main time "a_t": -T2 <= t <= 0; enforce
    // those constraints:
    assert(!IsPos(a_t));
    Time     t   = std::max(a_t, -m_T2);
    MassRate res = std::max(m_burnRateI2 + (m_bMu2 + m_aMu2 * t) * t,
                            MassRate(0.0));
    assert(IsFinite(res) && !IsNeg(res));   // MuHat1, must be > 0
    return res;
  }
  default:
    return MassRate(0.0);
  }
}

//===========================================================================//
// Angle-of-Attack:                                                          //
//===========================================================================//
Angle Ascent2::AoA(Time a_t) const
{
  switch (m_mode)
  {
  case FlightMode::Burn1:
  {
    // Arg: "tau" which is the time since Stage1 ignition (where the latter
    // is calculated via "m_T1", because the actual event-based "m_ignTime1"
    // is not known yet), so AoA1(tau=0)=0, ie @ launch;   0 <= tau <= T1;
    // NB: The coeff "b" @ "tau" has INVERTED sign:
    Time   tIgn1 = m_cutOffTime1  - m_T1;
    Time   tau   = std::min(std::max(a_t - tIgn1,    0.0_sec), m_T1);
    Angle  res   = std::max(tau * (m_aAoA1 * tau - m_bAoA1), 0.0_rad);
    if (res > 1.01 * MaxAoA1)
    {
      OutputCtls();
      throw std::logic_error
            ("Burn1: MaxAoA1=" + std::to_string(MaxAoA1.Magnitude()) +
             " exceeded: AoA=" + std::to_string(res    .Magnitude()));
    }
    return res;
  }

  case FlightMode::Gap:
    return 0.0_rad;

  case FlightMode::Burn2:
  {
    // Arg: "t" is time from Stage2 cut-off (same as orbital insertion
    // time), so it is just the main time "a_t": -T2 <= t <= 0; enforce
    // those constraints:
    assert(!IsPos(a_t));
    Time  t   = std::max(a_t, -m_T2);
    Angle res = std::max(t *  (m_aAoA2 * t + m_bAoA2), 0.0_rad);
    if (res > 1.01 * MaxAoA2)
    {
      OutputCtls();
      throw std::logic_error
            ("Burn2: MaxAoA2=" + std::to_string(MaxAoA2.Magnitude()) +
             " exceeded: AoA=" + std::to_string(res    .Magnitude()));
    }
    return res;
  }

  default:
    return 0.0_rad;
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
      Time     Isp = x * IspSL1 + (1.0 - x) * IspVac1;
      return a_burn_rate * Isp  * g0K;
    }

    case FlightMode::Gap:
      assert(IsZero(a_burn_rate));
      return ForceK(0.0);

    case FlightMode::Burn2:
      // We assume the Vacuum mode for Stage2:
      assert(IsPos(a_burn_rate));
      return a_burn_rate * IspVac2 * g0K;

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

  // Initialise "m" to the mass @ Orbital Insertion (incl the unspent
  // Stage2 Propellant):
  Mass m = m_emptyMass2 + m_propMass2 * PropRem2 + m_payLoadMass;

  // In any case, ADD the SpentPropMass (between "a_t" and Orbital Insert-
  // ion):
  m += spentPropMass;

  // If Stage1 has NOT separated yet, add its Empty and PropRemnants:
  if (m_mode != FlightMode::Burn2 && m_mode != FlightMode::Gap)
    m += m_emptyMass1 + m_propMass1 * PropRem1;

  // Add the FairingMass if the Fairing has NOT separated @  "a_t". In the
  // Bwd time, it means that "m_fairingSepTime" is known and "a_t" is below
  // it:
  if (!IsZero(m_fairingSepTime) && a_t < m_fairingSepTime)
    m += FairingMass;

  // All Done:
  assert(IsPos(m));
  return m;
}

//===========================================================================//
// "OutputCtls":                                                             //
//===========================================================================//
// For Testing Only:                                                       //
//
void Ascent2::OutputCtls() const
{
  if (m_os == nullptr)
    return;

  // Here "t" is the time from Stage2 cut-off (so t=-T2..0):
  std::cout << "# T2   := " << m_T2.Magnitude() << ';'    << std::endl;
  std::cout << "# AoA2 := t * ("   << m_aAoA2.Magnitude() << " * t + ("
            << m_bAoA2.Magnitude() << ")); "              << std::endl;
  std::cout << "# mu2  := " << m_burnRateI2.Magnitude()   << " + ("
            << m_bMu2.Magnitude()  << ") * t + ("
            << m_aMu2.Magnitude()  << ") * t^2;"          << std::endl;

  // Here "t" is the time from Stage1 cut-off (so t=-T1..0), and "tau" is
  // the time since Stage1 ignition (NB: the (-b) coeff at "tau" for AoA),
  // so tau=0..T1:
  std::cout << "# T1   := " << m_T1.Magnitude() << ';'    << std::endl;
  std::cout << "# AoA1 := tau * (" << m_aAoA1.Magnitude() << " * tau - ("
            << m_bAoA1.Magnitude() << ")); "              << std::endl;
  std::cout << "# mu1  := " << m_burnRateI1.Magnitude()   << " + ("
            << m_bMu1.Magnitude()  << ") * t + ("
            << m_aMu1.Magnitude()  << ") * t^2;"          << std::endl;
}

//===========================================================================//
// "GetRunRes":                                                              //
//===========================================================================//
// Extract the "RunRes" from the Cache or compute a new one. If either method
// fails, return the "nullopt":
//
std::option<Ascent2::RunRes> Ascent2::GetRunRes
(
  double const a_xs[NP],
  void*        a_env
)
{
  assert(a_env != nullptr);

  //-------------------------------------------------------------------------//
  // Construct the "AscCtls" obj used for caching:                           //
  //-------------------------------------------------------------------------//
  // a_xs = [aHat2, aMuHat2, aAoA2, bAoA2, TGap, aHat1, aMuHat1, aAoA1, bAoA1]:
  // XXX: Can we simply cast "a_xs" to an "AscCtls" ptr, w/o copying the data
  // over?
  AscCtls ctls
  {
    .m_aHat2  = a_xs[0], .m_muHat2 = a_xs[1],
    .m_aAoA2  = a_xs[2], .m_bAoA2  = a_xs[3],
    .m_TGap   = a_xs[4],
    .m_aHat1  = a_xs[5], .m_muHat1 = a_xs[6],
    .m_aAoA1  = a_xs[7], .m_bAoA1  = a_xs[8]
  };

  // Is the result already in the Cache?
  // # pragma omp critical(Cache)
  {
    auto it = s_Cache.find(ctls);
    if (it != s_Cache.end())
      return  it->second;
  }

  //-------------------------------------------------------------------------//
  // Otherwise, need a new Run:                                              //
  //-------------------------------------------------------------------------//
  // For safety, construct a new "Ascent2" obj from "a_env":
  Ascent2 const* proto = reinterpret_cast<Ascent2 const*>(a_env);
  Ascent2        asc(*proto);

  // Set the Ctl Params in "asc":
  asc.SetCtlParams(ctls);

  // Run the integrator:
  RunRes res  = asc.Run(); // NB: Exceptions are handled inside "Run"

  if (res.m_rc != RunRC::Error)
  {
    // Cache the results:
    // # pragma omp critical(Cache)
    s_Cache[ctls] = res;
    return res;
  }

  // If we got here: Both the Cache look-up and the new evaluation have
  // failed:
  return std::nullopt;
}

//===========================================================================//
// "EvalNLOptObjective":                                                     //
//===========================================================================//
// Evaluate the Function to be Minimised (by NLopt): It is actually the LV Mass
// at the end of Bwd integration:
//
double Ascent2::EvalNLoptObjective
(
  unsigned     DEBUG_ONLY(a_n),
  double const a_xs[NP],
  double*      DEBUG_ONLY(a_grad),
  void*        a_env
)
{
  assert(int(a_n) == NP && a_grad == nullptr && a_xs != nullptr);

  // Find or compute the "RunRes":
  std::optional<RunRes> res = GetRunRes(a_xs, a_env);

  return
    bool(res)
    ? // Got a valid result, use its mass:
      (res.value().m_mT)  .Magnitude()
    : // Otherwise, we return large invalid mass, but not too large (so not to
      // disrupt the optimisation algorithm):
      (2.0 * MaxStartMass).Magnitude();
}

//===========================================================================//
// "EvalNLOptConstraints":                                                   //
//===========================================================================//
// Evaluate the Constraints vector (all components must be <= 0) in NLopt:
//
void Ascent2::EvalNLOptConstraints
(
  unsigned     DEBUG_ONLY(a_m),
  double       a_constrs[],
  unsigned     DEBUG_ONLY(a_n),
  double const a_xs[NP],
  double*      DEBUG_ONLY(a_grad),
  void*        a_env
)
{
  assert(a_m == 2 && int(a_n) == NP && a_grad == nullptr && a_xs != nullptr);

  // Find or compute the "RunRes":
  std::optional<RunRes> res = GetRunRes(a_xs, a_env);

  if (bool(res))
  {
    a_constrs[0] = (std::max(res.value().m_hT, 0.0_km)   ).Magnitude();
    a_constrs[1] = (std::max(res.value().m_VT, VelK(0.0))).Magnitude();
  }
  else
  {
    // Failed to compute the "RunRes". Return positive but not too large
    // values:
    a_constrs[0] = (10.0_km).Magnitude();
    a_constrs[1] = VelK(1.0).Magnitude();
  }
}

//===========================================================================//
// "FindOptimalAscentCtls"                                                   //
//===========================================================================//
std::pair<Ascent2::AsctCtls, Mass> Ascent2::FindOptimalAscentCtls
(
  // Mission Params:
  Mass            a_payload_mass,
  LenK            a_h_perigee,
  LenK            a_h_apogee,
  Angle_deg       a_incl,
  Angle_deg       a_launch_lat,
  // LV Params (those which are considered to be non-"constexpr"):
  double          a_alpha1,      // FullMass1 / FullMass2
  ForceK          a_thrust2_vac,
  ForceK          a_thrust1_vac,
  // Logging Params:
  std::ostream*   a_os,
  int             a_log_level
)
{
  // Create the NLOpt obj. There are currently 9 params for optimization. Will
  // use the "COBYLA" method which can handle non-linear constraints  (for the
  // final "V" and "h" vals):
  nlopt::opt opt(nlopt::LN_COBYLA, NP);
  opt.set_min_objective(EvalNLOpt, this);

}

}
// End namespace SpaceBallistics
