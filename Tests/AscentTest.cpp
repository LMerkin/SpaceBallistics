// vim:ts=2:et
//===========================================================================//
//                          "Tests/AscentTest.cpp":                          //
//         Ascent to Low Earth Orbit Integration in Various Modes            //
//===========================================================================//
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include "SpaceBallistics/PhysEffects/LVAeroDyn.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"
#include <ostream>
#include <optional>
#include <utility>

using namespace  SpaceBallistics;
namespace EAM  = EarthAtmosphereModel;
namespace LVAD = LVAeroDyn;

namespace
{
  //=========================================================================//
  // The "Ascent" Class:                                                     //
  //=========================================================================//
  class Ascent
  {
  public:
    //-----------------------------------------------------------------------//
    // Consts:                                                               //
    //-----------------------------------------------------------------------//
    // Earth Params (TODO: Take the Flattening and J2 into account):
    constexpr static GMK      K              = BodyData<Body::Earth>::K;
    constexpr static AccK     g0K            = To_Len_km(g0);
    constexpr static LenK     R              = 6371.0_km;   // Equi-Volume
    // Launch Pad: Assume Lat=63 deg, and same Inclination.
    // Then the useful Earth Rotation Velocity is:
    constexpr static VelK     ERV            =
      TwoPi<double> * R / BodyData<Body::Earth>::SiderealRotationPeriod *
      Cos(To_Angle_rad(63.0_deg));

    // LV Params:
    // Stage1 (engine performance @ SL and in Vac); XXX: "BurnT1" may be for
    // reference only (since the mass rate may be variable); "PropMass1" is
    // quasi-fixed and plays the role of a scaling param of the whole stack:
    constexpr static Mass   PropMass1        = 228100.0_kg;
    constexpr static double K1               = 0.875;
    constexpr static Mass   EmptyMass1       = PropMass1 * (1.0 / K1 - 1.0);
    constexpr static Mass   FullMass1        = PropMass1 / K1;
    constexpr static Time   BurnT1           = 151.6_sec;
    constexpr static Time   IspSL1           = 326.0_sec;
    constexpr static Time   IspVac1          = 353.0_sec;
    constexpr static double PropReserve1     = 0.01;

    // Stage2 (only in Vac): Most params are dynamic:
    constexpr static double K2               = 0.9;
    constexpr static Time   IspVac2          = 381.0_sec;
    constexpr static double PropReserve2     = 0.01;

    // Total Start Mass Limit:
    constexpr static Mass   StartMassLimit   = 360000.0_kg;

    // Body and Fairing:
    // XXX: It is currently assumed that Fairing separation occurs whan the
    // dynamic pressure reaches 1 Pa or less:
    constexpr static Mass     FairingMass    = 1450.0_kg; // Soyuz-2.1b: 1100
    constexpr static Pressure FairingSepCond = Pressure(1.0);
    // LV Diameter:
    constexpr static LenK     D              = 0.0041_km;
    constexpr static decltype(Sqr(D)) CroS   = 0.25 * Pi<double> * Sqr(D);

    // ODE Integration Params: 1 msec step (can grow to 10 msec);
    // MaxStep correponds to roughly  80 m distance at the orbital velocity;
    // the relative error of 1e-6 correponds to ~7m distance and ~0.01 m/sec
    // velocity uncertainly at orbit:
    constexpr static Time   ODEInitStep      = 0.001_sec;
    constexpr static Time   ODEMaxStep       = 10.0 * ODEInitStep;
    constexpr static double ODERelPrec       = 1e-6;
    // Singular point detection:
    constexpr static VelK   SingV            = VelK(0.0001);    // 0.1 m/sec

    //-----------------------------------------------------------------------//
    // Types:                                                                //
    //-----------------------------------------------------------------------//
    // State Vector (for some time t <= 0, where t=t0=0 corresponds to the Orbi-
    // tal Insertion):
    // Computations are performed in TopoCentric (Start) Planar Polar CoOrds, so
    // the first 3 components are r, rDot, omega = phiDot. The last component is
    // the Total Spent Propellant Mass between the curr "t" and "t0" (which is a
    // continuous and DECREASING function of "t",  as opposed to the Total Mass
    // which is discontinuous when the Stage1 or Fairing are jettisones):
    //
    using StateV  = std::tuple<LenK, VelK, AngVel, Mass>;

    // The Time Derivative of the "StateV". The last component is the BurnRate:
    using DStateV = std::tuple<VelK, AccK, AngAcc, MassRate>;

    // Flight Mode:
    enum class FlightMode: int
    {
      UNDEFINED = 0,
      Burn1     = 1,  // Stage1 Burn,
      Gap       = 2,  // Ballistic Gap
      Burn2     = 3   // Stage2 Burn
    };

    inline char const* ToString(FlightMode a_mode)
    {
      switch (a_mode)
      {
        case FlightMode::Burn1: return "Burn1";
        case FlightMode::Gap  : return "Gap";
        case FlightMode::Burn2: return "Burn2";
        default:                return "UNDEFINED";
      }
    }

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // Performance Characteristics:
    LenK   const          m_hC;
    Mass   const          m_payLoadMass;

    // Stage2 Params:
    Mass   const          m_emptyMass2;
    Mass   const          m_propMass2;
    Mass   const          m_fullMass2;
    ForceK const          m_thrust2;

    // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
    // ignition):
    Time const            m_gapT;

    // Transient data (during flight path integration):
    FlightMode            m_mode;
    Time                  m_ignTime2;       // Burn2 start: St2 Ignition Time
    Time                  m_fairingSepTime; // Fairing Separation Time
    Time                  m_cutOffTime1;    // Gap   start: St1 Cut-Off  Time
    Time                  m_ignTime1;       // Burn1 start: St1 Ignition Time

    // Singular Point (if and when reached):
    std::optional<StateV> m_singS;
    std::optional<Time>   m_singT;

    // For output:
    std::ostream*         m_os;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // NB: We first detemined EmptyMass2 using the Nominal FullMass2 (w/o the
    // "a_prop_mass_adj2" coeff),  and then the Nominal and Actual PropMass2,
    // and then the Actual FullMass2; Thrust2 is determined using the TMR2 and
    // the Nominal FullMass2:
    Ascent
    (
      LenK          a_hc,
      Mass          a_payload_mass,
      double        a_mass12_ratio,
      double        a_tmr2,
      double        a_prop_mass_adj2,
      Time          a_gap_t,
      std::ostream* a_os
    )
    : m_hC            (a_hc),
      m_payLoadMass   (a_payload_mass),
      m_emptyMass2    (FullMass1    / a_mass12_ratio    * (1.0 - K2)),
      m_propMass2     (m_emptyMass2 / (1.0 / K2  - 1.0) * a_prop_mass_adj2),
      m_fullMass2     (m_emptyMass2 + m_propMass2),
      m_thrust2       (FullMass1    / a_mass12_ratio * a_tmr2 * g0K),
      m_gapT          (a_gap_t),
      m_mode          (FlightMode::UNDEFINED),
      m_ignTime2      (),
      m_fairingSepTime(),
      m_cutOffTime1   (),
      m_ignTime1      (),
      m_singS         (std::nullopt),
      m_singT         (std::nullopt),
      m_os            (a_os) // May be NULL
    {
      if (m_hC <=100.0_km     || IsNeg(m_payLoadMass) ||  IsNeg(m_gapT) ||
         !IsPos(m_emptyMass2) || IsNeg(m_propMass2)   || !IsPos(m_thrust2))
        throw std::invalid_argument("Ascent::Ctor: Invalid param(s)");

      // Check the MaxStartMass (the actual Start Mass may be lower if we find
      // that a lesser mass of Stage1 Propellant  is used up to the singular
      // point):
      Mass maxStartMass =
        FullMass1 + m_fullMass2 + m_payLoadMass + FairingMass;

      if (maxStartMass > StartMassLimit)
        throw std::invalid_argument("Ascent::Ctor: StartMassLimit Exceeded");
    }

    //-----------------------------------------------------------------------//
    // Integrate the Ascent Trajectory:                                      //
    //-----------------------------------------------------------------------//
    // Returns the altitude of the singular point (if reached), or "nullopt":
    //
    std::optional<std::pair<StateV, Time>> Run()
    {
      // We run the integration BACKWARDS from the orbital insertion point
      // (@ t=0):
      // Angular velocity at the Circular Orbit, XXX: with correction  for the
      // Earth Rotation Velocity (assuming the Target Orbit Inclination is the
      // same as the Launch Site Latitude):
      LenK   aC     = R + m_hC;
      VelK   v0     = SqRt(K / aC) - ERV;
      AngVel omega0 = 1.0_rad * v0 / aC;

      // The initial Radial Velocity is 0, and the final delta of Spent Mass
      // is also 0:
      StateV s0 = std::make_tuple(aC, VelK(0.0), omega0, 0.0_kg);

      // The RHS and the Call-Back Lambdas:
      auto rhs =
        [this](StateV const&  a_s, Time a_t) -> DStateV
        { return this->ODERHS(a_s, a_t); };

      auto cb  =
        [this](StateV* a_s,   Time a_t) -> bool
        { return this->ODECB (a_s, a_t); };

      // FlightMode ctl (modes are switched based on the SpentPropMass):
      m_mode           = FlightMode::Burn2;
      m_ignTime2       = 0.0_sec;  // Not known yet
      m_cutOffTime1    = 0.0_sec;  //
      m_fairingSepTime = 0.0_sec;  //
      m_ignTime1       = 0.0_sec;  //
      m_singS          = std::nullopt;
      m_singT          = std::nullopt;

      // Run the ODE Integrator. The maximum duration (which is certainly
      // enough for ascent to orbit) is 1 hour:
      constexpr Time t0    = 0.0_sec;
      constexpr Time tMin  = -3600.0_sec;

      // NB: All necessary exception handling is provided inside "RKF5":
      DEBUG_ONLY(Time tEnd =)
        RKF5(&s0, t0, tMin, rhs,
             -ODEInitStep, -ODEMaxStep, ODERelPrec, &cb, m_os);
      assert(tEnd >= tMin);

      // Have we reached a singular point (this is a necessary condition for a
      // successfult ascent to orbit)?
      assert(bool(m_singS) == bool(m_singT));
      return
        bool(m_singS)
        ? std::make_optional(std::make_pair(m_singS.value(), m_singT.value()))
        : std::nullopt;
    }

  private:
    //-----------------------------------------------------------------------//
    // ODE RHS for a simple Gravity Turn:                                    //
    //-----------------------------------------------------------------------//
    DStateV ODERHS(StateV const& a_s, Time a_t)
    {
      // NB: In the RHS evaluation, r <= R is allowed, as it does not cause any
      // singularities by itself; but we detect this condition in the Call-Back,
      // which means that the integration is over (successfully or otherwise):
      LenK     r       = std::get<0>(a_s);
      VelK     rDot    = std::get<1>(a_s);
      AngVel   omega   = std::get<2>(a_s);
      assert(!IsPos(a_t) && IsPos(r));

      // NB: omega < 0  is qualitatively impossible, because omega=0 implies
      // omegaDot  = 0; but it may occur due to a finite integration step,
      // so we have to control it manually:
      if (IsNeg(omega))
        omega = AngVel(0.0);

      // Atmospheric Conditions:
      // Curr altitude (possible slightly negative vals are rounded to 0):
      LenK   h   = std::max(r - R, 0.0_km);
      auto   atm = EAM::AirParams(h);

      // Propellang Burn Rate (>= 0):
      MassRate burnRate = PropBurnRate(a_t);
      assert(!IsNeg(burnRate));

      // Thrust:
      ForceK F  = Thrust(burnRate, std::get<0>(atm));
      assert(!IsNeg(F) && IsZero(F) == IsZero(burnRate));

      // Abs Velocity:
      VelK   V  =
        IsZero(omega)
        ? Abs (rDot)
        : SqRt(Sqr(rDot) + Sqr(r * omega / 1.0_rad));
      assert(!IsNeg(V));

      // The Mach Number:
      double M  = double(To_Len_m(V) / std::get<3>(atm));
      assert(M >= 0.0);

      // The Aerodynamic Force (assuming AngleOfAttack = 0):
      ForceK Q  =
        LVAD::cD(M, 0.0) * 0.5 * To_Len_km(std::get<1>(atm)) *
        Sqr(V)  * CroS;
      assert(!IsNeg(Q));

      Mass   m  = LVMass(a_s, a_t);
      assert(IsPos(m));

      // Acceleration (in the direction of the Velocity vector)
      // due to the Trust and AeroDynamic Force:
      AccK fqAcc = (F - Q) / m;

      // The RHS components:
      AccK   r2Dot;
      AngAcc omegaDot;

      if (LIKELY(!IsZero(fqAcc)))
      {
        // In this generic case, we need a non-0 V, otherwise the direction of
        // the vector "fqAcc" is not defined. Check that "V" is not too small:
        if (V < SingV)
        {
          m_singS =  a_s;
          m_singT =  a_t;
          throw std::logic_error("Singularity: V ~= 0");
        }

        r2Dot =
            r * Sqr(omega / 1.0_rad)   // "Kinematic"   term
          - K / Sqr(r)                 // Gravitational term
          + fqAcc * rDot  / V;         // "F" and "Q"   term

        omegaDot =
          omega *
          (
            - 2.0 * rDot / r           // "Kinematic"   term
            + fqAcc      / V           // "F" and "Q"   term
          );
      }
      else
      {
        // Degenerate case: no "F" and "Q" effects:
        r2Dot    = r * Sqr(omega / 1.0_rad) - K / Sqr(r);
        omegaDot = -2.0  * omega * rDot / r;
      }

      // The result: NB: the last derivative is <= 0, since the corresp
      // "StateV" components is int_t^0 burnRate(t') dt', ie "t" is the
      // LOWER integration limit:
      return std::make_tuple(rDot, r2Dot, omegaDot, -burnRate);
    }

    //-----------------------------------------------------------------------//
    // ODE CallBack (invoked after the completion of each RKF5 step):        //
    //-----------------------------------------------------------------------//
    // Here FlightMode switching occurs, so this method is non-"const":
    //
    bool ODECB(StateV* a_s, Time a_t)
    {
      assert(a_s != nullptr && !IsPos(a_t));
      LenK   r             = std::get<0>(*a_s);
      LenK   h             = r - R;
      VelK   rDot          = std::get<1>(*a_s);
      AngVel omega         = std::get<2>(*a_s);
      Mass   spentPropMass = std::get<3>(*a_s);
      Mass   m             = LVMass(*a_s, a_t); // Using the "old" Mode yet!

      if (IsNeg(omega))
      {
        // This is mathematically impossible (because omega=0 => dot_omega=0),
        // but may happen due to a finite time step. In that case, "manually"
        // correct it -- in the "StateV" as well:
        omega = AngVel(0.0);
        std::get<2>(*a_s) = omega;
      }

      // NB: "spentPropMass" devreases in time, so increases in Bwd time. It
      // is 0 at Orbit Isertion Time (t=0):
      assert(IsPos(r) && !IsNeg(spentPropMass) && IsPos(m));

      // Switching Burn2 -> Gap occurs according to "spentPropMass": when all
      // Stage2 propellants (except the Reserve) is spent:
      if (m_mode == FlightMode::Burn2  &&
          spentPropMass >= m_propMass2 * (1.0 - PropReserve2))
      {
        m_ignTime2    = a_t;
        m_mode        = FlightMode::Gap;
        // Furthermore, since the duration of the Ballistic Gap is known, at
        // this point we already know the Stage1 Cut-Off time:
        m_cutOffTime1 = a_t - m_gapT;
        if (m_os != nullptr)
          *m_os << "# t=" << a_t << ", H=" << h
                << ": Stage2 Ignition, Ballistic Gap Ends" << std::endl;
      }

      // Switching Gap -> Burn1 occurs merely by time delay (if the Gap is 0,
      // this will happen immediately after the above change Burn2 -> Gap):
      if (m_mode  == FlightMode::Gap)
      {
        assert(!IsZero(m_cutOffTime1));
        if (a_t <= m_cutOffTime1)
        {
          m_mode   = FlightMode::Burn1;
          if (m_os != nullptr)
            *m_os << "# t=" << a_t << ", H=" << h
                  << ": Stage1 Cut-Off, Ballistic Gap Starts" << std::endl;
        }
      }

      // Switching Burn1 -> UNDEFINED:
      // When we are in Stage1 burn and the whole expendable Propellant amt
      // has been spent:
      if (m_mode == FlightMode::Burn1 &&
          spentPropMass >= (PropMass1   * (1.0 - PropReserve1) +
                            m_propMass2 * (1.0 - PropReserve2)))
      {
        m_ignTime1 = a_t;
        m_mode     = FlightMode::UNDEFINED;
        if (m_os != nullptr)
          *m_os << "# t=" << a_t << ", H=" << h
                << ": Stage1 Ignition"   << std::endl;
      }

      // In any mode, detect the Fairing Separation Condition:
      // In the reverse time, it's when the Dynamic Pressure becomes HIGHER
      // that the threshold:
      auto     atm = EAM::AirParams(std::max(h, 0.0_km));
      Density  rho = std::get<1>(atm);
      auto     V2  = To_Len_m(Sqr(rDot) + Sqr(r * omega / 1.0_rad));
      Pressure Q   = 0.5 * rho * V2;

      if (IsZero(m_fairingSepTime) && Q >= FairingSepCond)
      {
        m_fairingSepTime = a_t;
        if (m_os != nullptr)
          *m_os << "# t=" << a_t << ", H=" << (r - R) << ": Fairing Separation"
                << std::endl;
      }
      // If we got into the UNDEFINED mode, stop now (all Propellant has been
      // exhaused, in the reverse time).  Also stop if we are on the surface:
      //
      bool cont = IsPos(h) && m_mode != FlightMode::UNDEFINED;

      // Abs Velocity: Similar to the "ODERHS",  we bound it away from 0  to
      // avoid singularities; "V" being close to 0 is NOT allowed because this
      // would make uncertain the direction of the Thrust vector (in the model
      // under consideration); V ~= 0 is formally acceptable in No-Thrust modes:
      if (cont && (m_mode == FlightMode::Burn1 || m_mode == FlightMode::Burn2))
      {
        VelK   V  =
          IsZero(omega)
          ? Abs (rDot)
          : SqRt(Sqr(rDot) + Sqr(r * omega / 1.0_rad));

        if (V < SingV) // As in the "ODERHS"
        {
          m_singS = *a_s;
          m_singT = a_t;
          cont    = false;
        }
      }

      // XXX: Output (with a 100 msec step, or if we are going to stop now):
      if (m_os  != nullptr &&
         (!cont || int(Round(double(a_t / 0.001_sec))) % 100 == 0))
      {
        // Thrust and BurnRate are for info only:
        MassRate burnRate = PropBurnRate(a_t);
        Pressure p        = std::get<0>(atm);
        auto     thrust   = Thrust(burnRate, p) / g0K;
        VelK     Vhor     = r * (omega / 1.0_rad);
        auto     V2       = Sqr(rDot) + Sqr(Vhor);
        Pressure q        = 0.5 * std::get<1>(atm) * To_Len_m(V2);
        VelK     V        = SqRt(V2);

        *m_os << a_t.Magnitude()      << '\t' << h.Magnitude()      << '\t'
              << rDot.Magnitude()     << '\t' << Vhor.Magnitude()   << '\t'
              << V.Magnitude()        << '\t' << m.Magnitude()      << '\t'
              << (m - m_payLoadMass).Magnitude()                    << '\t'
              << ToString(m_mode)     << '\t' << thrust.Magnitude() << '\t'
              << burnRate.Magnitude() << '\t' << q.Magnitude()
              << std::endl;
      }
      return cont;
    }

  public:
    //-----------------------------------------------------------------------//
    // Propellant Burn Rate:                                                 //
    //-----------------------------------------------------------------------//
    MassRate PropBurnRate(Time) const
    {
      // XXX: For the moment, the rates are constant, depend on the Mode only;
      // the exact Time is ignored. The value is >= 0:
      switch (m_mode)
      {
      case FlightMode::Burn1:
        // For Stage1, PropMass1 and BurnT1 are currently fixed:
        return PropMass1 / BurnT1;

      case FlightMode::Gap:
        return MassRate(0.0);

      case FlightMode::Burn2:
        // For Stage2, the BurnRate is also fixed, but via "m_thrust2":
        return m_thrust2 / (IspVac2 * g0K);

      default:
        return MassRate(0.0);
      }
    }

    //-----------------------------------------------------------------------//
    // Thrust:                                                               //
    //-----------------------------------------------------------------------//
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const
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
      __builtin_unreachable();
    }

    //-----------------------------------------------------------------------//
    // Current Mass:                                                         //
    //-----------------------------------------------------------------------//
    Mass LVMass(StateV const& a_s, Time a_t) const
    {
      // The mass of Propellants spent between "a_t" and the Orbital Insertion
      // instant (t=0):
      Mass spentPropMass = std::get<3>(a_s);
      assert(!IsNeg(spentPropMass));

      // Initialise "m" to the mass @ Orbital Insertion (incl the unspent
      // Stage2 Propellant):
      Mass m = m_emptyMass2 + m_propMass2 * PropReserve2 + m_payLoadMass;

      // In any case, ADD the SpentPropMass (between "a_t" and Orbital Insert-
      // ion):
      m += spentPropMass;

      // If Stage1 has NOT separated yet, add its Empty and PropReserve Masses:
      if (m_mode != FlightMode::Burn2 && m_mode != FlightMode::Gap)
        m += (EmptyMass1 + PropMass1 * PropReserve1);

      // Add the FairingMass if the Fairing has NOT separated @  "a_t". In the
      // Bwd time, it means that "m_fairingSepTime" is known and "a_t" is below
      // it:
      if (!IsZero(m_fairingSepTime) && a_t < m_fairingSepTime)
        m += FairingMass;

      // All Done:
      assert(IsPos(m));
      return m;
    }
  };

  //=========================================================================//
  // "FindLowestSingularPoint":                                              //
  //=========================================================================//
  // Find the maximum PayLoadMass such that the singular  point is reachable
  // (so the corresp Altitude is at minimum >= 0);
  // returns (PayLoadMass, SingularPointAltitude, StartMass). If the singular
  // point is not reachable, returns "nullopt":
  //
  std::tuple<Mass, LenK, Mass> FindLowestSingularPoint
  (
    LenK   a_hc,
    double a_mass12_ratio,
    double a_tmr2,
    double a_prop_mass_adj2,
    Time   a_gap_t
  )
  {
    // The result is accumulated here:
    Mass mL0        = 0.0_kg;
    LenK h0           (Inf<double>);
    Mass startMass0 = 0.0_kg;

    Mass mL         = mL0;         // Searched upwards
    Mass dmL        = 1000.0_kg;   // Later may be reduced

    while (true)
    {
      // Create an "Ascent" object and run the trajectory integration:
      Ascent asc(a_hc,    mL, a_mass12_ratio, a_tmr2, a_prop_mass_adj2,
                 a_gap_t, nullptr);

      auto   res = asc.Run();

      // If we have not reached the singular point, the upper bound has been
      // encountered:
      if (!bool(res))
      {
        // "mL" is beyond the valid range, so the currently-achieved "mL0" is
        // kept; continue searching upwards from it, with a reduced step.
        // However, we may also want to stop now:
        if (dmL <= 10.0_kg)
        {
          // We have already obtained "mL0" with sufficient accuracy; extract
          // the corresp "h" which must be >= 0  (otherwise, we would not get
          // to the singular point in the first place). If, however, "mL0" is
          // 0, the singular point is not reachable at all:
          if (UNLIKELY(IsZero(mL0)))
            throw std::logic_error("FindLowestSingularPoint: UnReachable");

          // If "mL0" is valid: All Done:
          return std::make_tuple(mL0, h0, startMass0);
        }
        // Otherwise, continue seraching from the unchanged "mL0" upwards,
        // with a reduced step:
        dmL /= 10.0;
        mL   = mL0 + dmL;
      }
      else
      {
        // Otherwise, "res" is valid, so update the curr result (after "Run"):
        mL0                      = mL;
        Ascent::StateV const& s0 = res.value().first;
        Time                  t0 = res.value().second;
        h0                       = std::get<0>(s0) - Ascent::R;
        assert(!IsNeg(h0));
        startMass0               = asc.LVMass(s0, t0);

        // Continue with the same step:
        mL               = mL0 + dmL;
      }
    }
    __builtin_unreachable();
  }

  //=========================================================================//
  // "FindOptTGap":                                                          //
  //=========================================================================//
  // The control parameter is "TGap". Find a suitable value of "TGap" such that
  // the singular point is at or below the 50 m elevation  (for technical reas-
  // ons, we cannot make it exactly 0). Returns (mL, TGap, StartMass):
  //
  std::tuple<Mass, Time, Mass> FindOptTGap
  (
    LenK   a_hc,
    double a_mass12_ratio,
    double a_tmr2,
    double a_prop_mass_adj2,
    bool   a_verbose
  )
  {
    if (a_verbose)
      std::cout << "hC=" << a_hc.Magnitude()      << ", Mass1/Mass2="
                << a_mass12_ratio    << ", TMR2=" << a_tmr2
                << ", PropMass2Adj=" << a_prop_mass_adj2 << std::endl;

    // XXX: Very primitive "brute-force" minimisation of "h0":
    // Typically, with "TGap" increasing, "h0" decreases at the rate of ~1 km
    // per 1 "TGap" sec, but "mL" also decreases; then a region of near-0 "h0"
    // is reached  (for a continuous range of "TGap"s).
    // So we stop at the LOWEST admissible "TGap". Initial estimate: TGap = 0:
    //
    auto res0 =
      FindLowestSingularPoint
        (a_hc,  a_mass12_ratio, a_tmr2, a_prop_mass_adj2, 0.0_sec);
    LenK h0   = std::get<1>(res0);

    // Then we can guess the initial value and the initial step of "TGap":
    Time tg0    = Round(0.5 * h0 / VelK(1.0));
    Time tgStep = (tg0 > 20.0_sec) ? 5.0_sec : 1.0_sec;

    // March with the initial "tgStep" until we get (h < 0.5_km); it's an error
    // if we could not get it:
    constexpr Time MaxTGap = 250.0_sec;
    for (Time tg = tg0; tg < MaxTGap; tg += tgStep)
    {
      auto res =
        FindLowestSingularPoint
          (a_hc, a_mass12_ratio, a_tmr2, a_prop_mass_adj2, tg);
      Mass mL        = std::get<0>(res);
      LenK h         = std::get<1>(res);
      Mass startMass = std::get<2>(res);

      if (a_verbose)
        std::cout << "\tTGap=" << tg.Magnitude() << ", h=" << h.Magnitude()
                  << ", mL="   << mL.Magnitude() << ", startMass="
                  << startMass.Magnitude()       << std::endl;

      // If we got VERY close to h=0, we are done:
      if (h < 0.075_km)
        return std::make_tuple(mL, tg, startMass);

      // Otherwise: If we still got sufficiently close to h=0, reduce the step
      // and proceed further:
      if (h < 0.5_km)
        tgStep = 0.1_sec;
    }
    // XXX: If we got here, an admissible singular point has not been found:
    throw std::logic_error("FindOptTGap: Not Found");
  }

  //=========================================================================//
  // "Stage2Options":                                                        //
  //=========================================================================//
  void Stage2Options(LenK a_hc)
  {
#   pragma omp parallel for collapse(3)
    for (int mass12Rpct = 300;   mass12Rpct <= 450;   mass12Rpct += 5)
    for (int tmr2pct    =  90;   tmr2pct    <= 190;   tmr2pct    += 5)
    for (int adj2pct    = 100;   adj2pct    <= 105;   adj2pct    += 1)
    {
      double mass12R    = double(mass12Rpct) / 100.0;
      double tmr2       = double(tmr2pct)    / 100.0;
      double adj2       = double(adj2pct)    / 100.0;
      try
      {
        auto [mL, TGap, startMass] =
          FindOptTGap(a_hc, mass12R, tmr2, adj2, false);

#       pragma omp critical
        std::cout
          << a_hc.Magnitude()      << '\t' << mass12R          << '\t'
          << tmr2                  << '\t' << adj2             << '\t'
          << TGap.Magnitude()      << '\t' << mL.Magnitude()   << '\t'
          << startMass.Magnitude() << '\t' << double(mL/startMass)
          << std::endl;
      }
      catch (...) {}
    }
  }
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main(int argc, char* argv[])
{
  using namespace std;

  if (argc < 2)
  {
    cerr << "PARAMETER: LEO_Altitude_km" << endl;
    return 1;
  }
  LenK          hC { atof(argv[1]) };
  Stage2Options(hC);
	return  0;
}
