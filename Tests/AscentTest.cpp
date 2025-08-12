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
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // Earth Params (TODO: Take the Flattening and J2 into account):
    constexpr static GMK      K              = BodyData<Body::Earth>::K;
    constexpr static LenK     R              = 6371.0_km;   // Equi-Volume
    // Launch Pad: Assume Lat=63 deg, and same Inclination.
    // Then the useful Earth Rotation Velocity is:
    constexpr static VelK     ERV            =
      TwoPi<double> * R / BodyData<Body::Earth>::SiderealRotationPeriod *
      Cos(To_Angle_rad(63.0_deg));

    // LV Params:
    // Stage1 (engine performance @ SL and in Vac); its full mass is determined
    // dynamically, but we need some initial estimates for reference:
    constexpr static double K1               = 0.9;
    constexpr static Time   IspSL1           = 326.0_sec;
    constexpr static Time   IspVac1          = 353.0_sec;
    constexpr static double PropRem1         = 0.01;

    // Stage2 (only in Vac): Most params are dynamic. Here "PropRem2" is a frac-
    // tion of the total Stage2 propellant load, not a fixed value:
    constexpr static double K2               = 0.93;
    constexpr static Time   IspVac2          = 381.0_sec;
    constexpr static double PropRem2         = 0.01;

    // Total Start Mass (Constant):
    constexpr static Mass   StartMass        = 360'000.0_kg;

    // Body and Fairing:
    // XXX: It is currently assumed that Fairing separation occurs whan the
    // dynamic pressure reaches 1 Pa or less:
    constexpr static Mass     FairingMass    = 1000.0_kg; // Soyuz-2.1b: 1100
    constexpr static Pressure FairingSepCond = Pressure(1.0);
    // LV Diameter:
    constexpr static LenK     D              = To_Len_km(4.1_m);
    constexpr static decltype(Sqr(D)) CroS   = 0.25 * Pi<double> * Sqr(D);

    // ODE Integration Params: 1 msec step; it may only be reduced, never
    // increased beyond the original value:
    constexpr static Time   ODEInitStep      = 0.001_sec;
    constexpr static Time   ODEMaxStep       = ODEInitStep;
    constexpr static double ODERelPrec       = 1e-6;

    // Singular point detection criteria: NB: "Vhor" approaches 0 much faster
    // than "Vr":
    constexpr static VelK   SingVr           = VelK(1e-2); // 10   m/sec
    constexpr static VelK   SingVhor         = VelK(1e-4); //  0.1 m/sec

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
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
    // Exception thrown when the flight path is approaching the singular point:
    class NearSingularityExn
    {
    public:
      // Data Flds:
      LenK  const m_r;
      VelK  const m_Vr;
      Mass  const m_spentPropMass;
      Time  const m_t;

      // Non-Default Ctor:
      NearSingularityExn(LenK a_r, VelK a_vr, Mass a_spent_prop_mass, Time a_t)
      : m_r             (a_r),
        m_Vr            (a_vr),
        m_spentPropMass (a_spent_prop_mass),
        m_t             (a_t)
      {}
    };

    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    Mass     const        m_payLoadMass;

    // Stage2 Params:
    Mass     const        m_fullMass2;
    Mass     const        m_emptyMass2;
    Mass     const        m_propMass2;
    ForceK   const        m_thrustVac2;
    MassRate const        m_burnRate2;

    // Stage1 Params:
    Mass     const        m_fullMass1;
    Mass     const        m_emptyMass1;
    Mass     const        m_propMass1;
    ForceK   const        m_thrustVac1;
    MassRate const        m_burnRate1;

    // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
    // ignition):
    Time                  m_TGap;

    // Transient data (during flight path integration): St2 Cut-Off Time is 0:
    FlightMode            m_mode;
    Time                  m_ignTime2;       // Burn2 start: St2 Ignition Time
    Time                  m_fairingSepTime; // Fairing Separation Time
    Time                  m_cutOffTime1;    // Gap   start: St1 Cut-Off  Time
    Time                  m_ignTime1;       // Burn1 start: St1 Ignition Time

    // Singular Point (if and when reached):
    std::optional<StateV> m_singS;
    std::optional<Time>   m_singT;

    // For output:
    std::ostream* const   m_os;
    int           const   m_logLevel;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    Ascent
    (
      Mass          a_payload_mass,
      double        a_alpha1,      // FullMass1 / FullMass2
      ForceK        a_thrust2_vac,
      ForceK        a_thrust1_vac,
      std::ostream* a_os,
      int           a_log_level
    )
    : m_payLoadMass   (a_payload_mass),
      m_fullMass2     ((StartMass - FairingMass - m_payLoadMass) /
                       (1.0 + a_alpha1)),
      m_emptyMass2    (m_fullMass2   * (1.0 - K2)),
      m_propMass2     (m_fullMass2   * K2),
      m_thrustVac2    (a_thrust2_vac),
      m_burnRate2     (m_thrustVac2  / (IspVac2  * g0K)),
      m_fullMass1     ((StartMass - FairingMass - m_payLoadMass) /
                       (1.0 + a_alpha1) * a_alpha1),
      m_emptyMass1    (m_fullMass1  * (1.0 - K1)),
      m_propMass1     (m_fullMass1  * K1),
      m_thrustVac1    (a_thrust1_vac),
      m_burnRate1     (m_thrustVac1  / (IspVac1  * g0K)),
      m_TGap          (),     // 0 by default
      m_mode          (FlightMode::UNDEFINED),
      m_ignTime2      (),
      m_fairingSepTime(),
      m_cutOffTime1   (),
      m_ignTime1      (),
      m_singS         (std::nullopt),
      m_singT         (std::nullopt),
      m_os            (a_os), // May be NULL
      m_logLevel      (a_log_level)
    {
      if ( IsNeg(m_payLoadMass) ||  IsNeg(m_TGap)      || !IsPos(m_fullMass2) ||
          !IsPos(m_thrustVac2)  || !IsPos(m_fullMass1) || !IsPos(m_thrustVac1))
        throw std::invalid_argument("Ascent::Ctor: Invalid param(s)");

      // Therefore:
      assert(IsPos(m_emptyMass2) && IsPos(m_propMass2) && IsPos(m_burnRate2) &&
             IsPos(m_emptyMass1) && IsPos(m_propMass1) && IsPos(m_burnRate1));
    }

    //=======================================================================//
    // "Run": Integrate the Ascent Trajectory:                               //
    //=======================================================================//
    // If the Singular Point has been reached, returns (FinalState, FinalTime);
    // otherwise, returns "nullopt":
    //
    std::optional<std::pair<StateV, Time>> Run(LenK a_hc, Time a_t_gap)
    {
      if (a_hc < 100.0_km || IsNeg(a_t_gap))
        throw std::invalid_argument("Assent::Run: Invalid hC or TGap");

      // We run the integration BACKWARDS from the orbital insertion point
      // (@ t=0):
      // Angular velocity at the Circular Orbit, XXX: with correction  for the
      // Earth Rotation Velocity (assuming the Target Orbit Inclination is the
      // same as the Launch Site Latitude):
      LenK   aC     = R + a_hc;
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
      m_TGap           = a_t_gap;
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
      try
      {
        DEBUG_ONLY(Time tEnd =)
          RKF5(&s0, t0, tMin, rhs,
               -ODEInitStep, -ODEMaxStep, ODERelPrec, &cb, m_os);
        assert(tEnd >=  tMin);
      }
      catch (NearSingularityExn const& ns)
      {
        // We have reached a vicinity of the singular point, this is a NORMAL
        // (and moreover,  a desirable) outcome. Compute a more precise sing-
        // ular point position:
        auto optSing = LocateSingularPoint(ns);
        if (bool(optSing))
        {
          m_singS = optSing.value().first;
          m_singT = optSing.value().second;
        }
      }
      catch (std::exception const& exn)
      {
        // Any other "standard" exceptions: Ascent unsuccessful:
        if (m_os != nullptr && m_logLevel >= 1)
          *m_os  << "# Ascent::Run: Exception: " << exn.what() << std::endl;
        return std::nullopt;
      }

      // Have we reached the singular point (this is a necessary condition for
      // a successfult ascent to orbit)?
      assert(bool(m_singS) == bool(m_singT));
      return
        bool(m_singS)
        ? std::make_optional(std::make_pair(m_singS.value(), m_singT.value()))
        : std::nullopt;
    }

  private:
    //=======================================================================//
    // ODE RHS for a simple Gravity Turn:                                    //
    //=======================================================================//
    DStateV ODERHS(StateV const& a_s, Time a_t)
    {
      // NB: In the RHS evaluation, r <= R is allowed, as it does not cause any
      // singularities by itself; but we detect this condition in the Call-Back,
      // which means that the integration is over (successfully or otherwise):
      LenK   r             = std::get<0>(a_s);
      VelK   Vr            = std::get<1>(a_s);
      AngVel omega         = std::get<2>(a_s);
      Mass   spentPropMass = std::get<3>(a_s);
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
        throw NearSingularityExn(r, Vr, spentPropMass, a_t);

      // Generic Case:
      // Abs Velocity: Since we are NOT near the singularity, it should be
      // sufficiently far away from 0:
      VelK   V  = SqRt(Sqr(Vr) + Sqr(Vhor));
      assert(IsPos(V));

      // AeroDynamic Drag:
      auto [atm, drag]  = AeroDynForce(r, V);
      assert(!IsNeg(drag));

      // Propellant Burn Rate (>= 0) and Thrust:
      MassRate burnRate = PropBurnRate(a_t);
      assert(!IsNeg(burnRate));

      ForceK   thrust   = Thrust(burnRate, std::get<0>(atm));
      assert(!IsNeg(thrust) && IsZero(thrust) == IsZero(burnRate));

      // The curr mass:
      // XXX: If it exceeds the over-all limit, there is no point in continuing
      // the Bwd intrgration; singular point is NOT reached in that case:
      Mass   m  = LVMass(a_s, a_t);
      assert(IsPos(m));

      // Acceleration (in the direction of the Velocity vector)
      // due to the Trust and AeroDynamic Drag Force:
      AccK tdAcc = (thrust - drag) / m;

      // The RHS components:
      AccK   r2Dot;
      AngAcc omegaDot;

      r2Dot =
          r * Sqr(omega / 1.0_rad)   // "Kinematic"     term
        - K / Sqr(r)                 // Gravitational   term
        + tdAcc * Vr    / V;         // Thrust and Drag term

      omegaDot =
        omega *
        (
          - 2.0 * Vr    / r          // "Kinematic"     term
          + tdAcc       / V          // Thrust and Drag term
        );

      // The result: NB: the last derivative is <= 0, since the corresp
      // "StateV" components is int_t^0 burnRate(t') dt', ie "t" is the
      // LOWER integration limit:
      return std::make_tuple(Vr, r2Dot, omegaDot, -burnRate);
    }

    //=======================================================================//
    // ODE CallBack (invoked after the completion of each RKF5 step):        //
    //=======================================================================//
    // Here FlightMode switching occurs, so this method is non-"const":
    //
    bool ODECB(StateV* a_s, Time a_t)
    {
      assert(a_s != nullptr && !IsPos(a_t));
      LenK   r             = std::get<0>(*a_s);
      LenK   h             = r - R;
      VelK   Vr            = std::get<1>(*a_s);
      AngVel omega         = std::get<2>(*a_s);
      Mass   spentPropMass = std::get<3>(*a_s);
      Mass   m             = LVMass(*a_s, a_t); // Using the "old" Mode yet!

      //---------------------------------------------------------------------//
      // Similar to "ODERHS", check if we are approaching the singularity:   //
      //---------------------------------------------------------------------//
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
        throw NearSingularityExn(r, Vr, spentPropMass, a_t);

      //---------------------------------------------------------------------//
      // Generic Case:                                                       //
      //---------------------------------------------------------------------//
      // NB: "spentPropMass" decreases over time, so increases in Bwd time.
      // It is 0 at Orbit Insertion Time (t=0):
      assert(IsPos(r) && !IsNeg(spentPropMass) && IsPos(m));

      // Output at the beginning:
      if (IsZero(a_t) && m_os != nullptr && m_logLevel >= 1)
      {
        assert(m_mode == FlightMode::Burn2);
        *m_os << "# t="    << a_t.Magnitude() << " sec, H=" << h.Magnitude()
              << " km, V=" << V.Magnitude()   << " km/sec, Mass="
              << m.Magnitude()                << " kg"      << std::endl;
      }
      //---------------------------------------------------------------------//
      // Switching Burn2 -> Gap:                                             //
      //---------------------------------------------------------------------//
      // Occurs according to "spentPropMass": when all Stage2 propellants
      // (except the Remnants) is spent:
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
        assert(m_curOffTime1 <= m_ignTime2);

        if (m_os != nullptr && m_logLevel >= 1)
          *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
                << " km, V=" << V.Magnitude()   << " km/sec, m="
                << m.Magnitude()                << " kg: "
                   "Stage2 Ignition, Ballistic Gap Ends"      << std::endl;
      }
      //---------------------------------------------------------------------//
      // Switching Gap -> Burn1:                                             //
      //---------------------------------------------------------------------//
      // Occurs merely by time delay (if the Gap is 0, this will happen immedi-
      // ately after the above change Burn2 -> Gap):
      //
      if (m_mode  == FlightMode::Gap)
      {
        assert(IsNeg(m_cutOffTime1));
        if (a_t <= m_cutOffTime1)
        {
          m_mode   = FlightMode::Burn1;

          // Re-calculate the Mass in Burn1 for the output (with Stage1):
          m        = LVMass(*a_s, a_t);
          if (m_os != nullptr && m_logLevel >= 1)
            *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
                  << " km, V=" << V.Magnitude()   << " km/sec, m="
                  << m.Magnitude()                << " kg: "
                     "Stage1 Cut-Off,  Ballistic Gap Starts"    << std::endl;
        }
      }
      //---------------------------------------------------------------------//
      // Switching Burn1 -> UNDEFINED:                                       //
      //---------------------------------------------------------------------//
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
          *m_os << "# t="    << a_t.Magnitude() << " sec, h=" << h.Magnitude()
                << " km, V=" << V.Magnitude()   << " km/sec, m="
                << m.Magnitude()                << " kg: "
                   "Stage1 Ignition"            << std::endl;
      }

      //---------------------------------------------------------------------//
      // In any mode, detect the Fairing Separation Condition:               //
      //---------------------------------------------------------------------//
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
          *m_os << "# t=" << a_t.Magnitude() << " sec, h="
                << (r - R).Magnitude()       << " km: Fairing Separation"
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
          *m_os << "# t=" << a_t.Magnitude()            << " sec: H=0, Vr="
                << Vr.Magnitude()   << " km/sec, Vhor=" << Vhor.Magnitude()
                << " km/sec, m="    << m.Magnitude()    << " kg, "
                << ToString(m_mode) << std::endl;
      }
      if (m_mode == FlightMode::UNDEFINED)
      {
        cont = false;

        if (m_os != nullptr && m_logLevel >= 1)
          *m_os << "# t=" << a_t.Magnitude()            << " sec: UNDEF, Vr="
                << Vr.Magnitude()   << " km/sec, Vhor=" << Vhor.Magnitude()
                << " km/sec, m="    << m.Magnitude()    << " kg, "
                << ToString(m_mode) << std::endl;
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

    //=======================================================================//
    // Atmospheric Conditions and Aerodynamic Drag Force:                    //
    //=======================================================================//
    static std::pair<EAM::AtmConds, ForceK>AeroDynForce(LenK a_r, VelK a_v)
    {
      // Curr altitude (possible slightly negative vals are rounded to 0):
      auto atm  = EAM::GetAtmConds(a_r - R);

      // The Mach Number:
      double M  = double(To_Len_m(a_v) / std::get<3>(atm));
      assert(M >= 0.0);

      // The Aerodynamic Drag Force (assuming AngleOfAttack = 0):
      ForceK drag =
        LVAD::cD(M, 0.0) * 0.5 * To_Len_km(std::get<1>(atm)) *
        Sqr(a_v) * CroS;

      return std::make_pair(atm, drag);
    }

    //=======================================================================//
    // "LocateSingularPoint":                                                //
    //=======================================================================//
    // The final stage of integration, where we assume omega=0 (a purely vert-
    // ical motion)   and integrate the simplified ODEs ANALYTICALLY to avoid
    // numerical instabilities in the viciniy of the singular point.
    // Returns (singH, singT)  if the singular point has been found, otherwise
    // "nullopt":
    //
    std::optional<std::pair<StateV, Time>>
    LocateSingularPoint(NearSingularityExn const& a_nse)
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
        *m_os << "Ascent::LocateSingularPoint: WARNING: h="
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

      // The Curr LV Mass:
      StateV s1    { r1, Vr1, AngVel(0.0), a_nse.m_spentPropMass };
      Mass   m1  = LVMass(s1, t1);
      assert(IsPos(m1));

      // The (constant) acceleration we will use:
      AccK   acc1 = thrust1 / m1 - g1;
      if (!IsPos(acc1))
      {
        // Then the LV is falling back to the pad, and the singular point is
        // not reachable; this is because we have arrived at the mass "m1"
        // which is too large:
        if (m_os != nullptr && m_logLevel >= 2)
          *m_os << "# Ascent::LocateSingularPoint: UnReachable: Mass=Acc="
                << (double(acc1 / g1) - 1.0) << " g" << std::endl;
        return std::nullopt;
      }

      // Otherwise: Remaining Bwd Time and Distance to the singular point:
      Time tau = Vr1            / acc1;
      LenK dr  = 0.5 * Sqr(Vr1) / acc1;
      assert(!(IsNeg(tau) || IsNeg(dr)));

      // Finally: State and Time of the singular point:  XXX: We do not check
      // for rS < R, it is acceptable (the neg value would be very small, any-
      // way):
      LenK rS    = r1 - dr;
      Time tS    = t1 - tau;
      Mass propS = a_nse.m_spentPropMass + burnRate1 * tau;
      StateV singS {rS, VelK(0.0), AngVel(0.0), propS};

      if (m_os != nullptr && m_logLevel >= 2)
        *m_os << "# SingularPoint Located: t1="    << t1.Magnitude()
              << " sec, tau="   << tau.Magnitude() << " sec, tS="
              << tS.Magnitude() << " sec, m1="     << m1.Magnitude()
              << " kg, mS="     << (m1 + burnRate1 * tau).Magnitude()
              << " kg, mS'="    << LVMass(singS, tS).Magnitude()
              << std::endl;

      return std::make_optional(std::make_pair(singS, tS));
    }

  public:
    //=======================================================================//
    // Propellant Burn Rate:                                                 //
    //=======================================================================//
    MassRate PropBurnRate(Time) const
    {
      // XXX: For the moment, the rates are constant, depend on the Mode only;
      // the exact Time is ignored. The value is >= 0:
      switch (m_mode)
      {
      case FlightMode::Burn1:
        return m_burnRate1;

      case FlightMode::Gap:
        return MassRate(0.0);

      case FlightMode::Burn2:
        return m_burnRate2;

      default:
        return MassRate(0.0);
      }
    }

    //=======================================================================//
    // Thrust:                                                               //
    //=======================================================================//
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
    }

    //=======================================================================//
    // Current Mass:                                                         //
    //=======================================================================//
    Mass LVMass(StateV const& a_s, Time a_t) const
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

    //=======================================================================//
    // "FindFeasibleTrajectory":                                             //
    //=======================================================================//
    void FindFeasibleTrajectory(LenK a_hc)
    {
      // Have we got at least one singular point?
      bool gotS = false;

      for (Time TGap = 0.0_sec; TGap <= 420.0_sec; TGap += 1.0_sec)
      {
        if (m_os != nullptr && m_logLevel >= 1)
          *m_os << "\n=====> TGap=" << TGap.Magnitude() << std::endl;

        // Run the integrator and hopr to reach the singular point:
        auto res = Run(a_hc, TGap);
        if (!res)
        {
          // If there was already a singular point, but not anymore, no point
          // in iterating further:
          if (gotS)
            break;
          else
            continue;  // Try, try again!
        }
        // If we got here: singular point reached:
        gotS    = true;
        LenK hS = std::get<0>(res.value().first) - R;

        if (m_os != nullptr)
          *m_os << "-----> hS=" << hS.Magnitude() << " km" << std::endl;
      }
    }
  };
  // End of "Ascent" Class

/*
  //=========================================================================//
  // "FindLowestSingularPoint":                                              //
  //=========================================================================//
  // Tries to find the maximum Stage2 Mass   (corresp to the minimim Singular
  // Point Altitude), if the signular point exists at all. In this case,  re-
  // turns (Stage2Mass, SingularPointAltitude, StartMass). Otherwise, returns
  // "nullopt":
  //
  void FindLowestSingularPoint
  (
    LenK   a_hc,
    Mass   a_payload_mass,
    ForceK a_thrust2_vac,
    ForceK a_thrust1_vac,
    Time   a_gap_t
  )
  {
    // The last "m2" mass for which the signular point has been reached:
    Mass m2SP = 0.0_kg;

    for (Mass m2 = 100'000.0_kg; m2 >= 20'000.0_kg; m2 -= 250.0_kg)
    {
std::cout << "==> m2 = " << m2.Magnitude() << " kg" << std::endl;
      // XXX: At this point, we don't have the Stage1 Mass, so use an estimate.
      // If the singular point has been reached, we do further iterations with
      // more precise Stage1 Mass -- and hope that the singular point will per-
      // sist:
      std::optional<std::pair<Ascent::StateV, Time>> res = std::nullopt;
      try
      {
        // Create an "Ascent" object and run the trajectory integration:
        Ascent asc(a_hc,    a_payload_mass, m2,      a_thrust2_vac,
                   a_thrust1_vac,  a_gap_t, &std::cout, 1);

        // NB: any exceptions apart from "NearSingularExn" are propagated to
        // here:
        res  = asc.Run();

        // So: have we got to the singular point?
        if (bool(res))
        {
          // Yes, singular point has been found:
          Mass startMass = asc.LVMass(res.value().first, res.value().second);
          Mass m1        = startMass - a_payload_mass - m2 -
                           Ascent::FairingMass;
          assert(IsPos(m1));

          m2SP                    = m2;
          Ascent::StateV const& s = res.value().first;
          LenK                  h = std::get<0>(s) - Ascent::R;
          assert(IsPos(m2SP) && !IsNeg(h));

Time ts = res.value().second;
std::cout << "SP: t = "       << ts.Magnitude()
          << " sec, m2 = "    << m2.Magnitude()
          << " kg, m1 = "     << m1.Magnitude()
          << " kg, m0 = "     << startMass.Magnitude()
          << " kg, h = "      << h.Magnitude()
          << " km"            << std::endl;
        }
        else
        if (!IsZero(m2SP))
          // If there was already a Singular Point encountered, but not with
          // this "m2", stop further searching -- we have left the feasibility
          // region:
          break;
      }
      catch (...)
      {
        // Any errors: "res" remains "nullopt", so this run is unsuccessful:
        assert(!bool(res));
      }
    }
    // End of "m2" loop
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
    double a_tmr2,
    bool   a_verbose
  )
  {
    if (a_verbose)
      std::cout << "hC=" << a_hc.Magnitude()      <<
                << ", TMR2=" << a_tmr2 << std::endl;

    // XXX: Very primitive "brute-force" minimisation of "h0":
    // Typically, with "TGap" increasing, "h0" decreases at the rate of ~1 km
    // per 1 "TGap" sec, but "mL" also decreases; then a region of near-0 "h0"
    // is reached  (for a continuous range of "TGap"s).
    // So we stop at the LOWEST admissible "TGap". Initial estimate: TGap = 0:
    //
    auto res0 =
      FindLowestSingularPoint(a_hc, a_mass12_ratio, a_tmr2, 0.0_sec);
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
        FindLowestSingularPoint(a_hc, a_mass12_ratio, a_tmr2, tg);
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
#   pragma omp parallel for collapse(2)
    for (int tmrSpct= 300;  mass12Rpct <= 450; mass12Rpct += 5)
    for (int tmr2pct=  90;  tmr2pct    <= 190; tmr2pct    += 5)
    {
      double mass12R    = double(mass12Rpct) / 100.0;
      double tmr2       = double(tmr2pct)    / 100.0;
      try
      {
        auto [mL, TGap, startMass] =
          FindOptTGap(a_hc, mass12R, tmr2, false);

#       pragma omp critical
        std::cout
          << a_hc.Magnitude()      << '\t' << mass12R          << '\t'
          << tmr2                  << '\t' << TGap.Magnitude() << '\t'
          << mL.Magnitude()        << '\t'
          << startMass.Magnitude() << '\t' << double(mL/startMass)
          << std::endl;
      }
      catch (...) {}
    }
  }
*/
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main(int argc, char* argv[])
{
  using namespace std;

  if (argc < 4)
  {
    cerr << "PARAMETERS: LEO_Altitude_km PayLoad_kg m1/m2_Ratio" << endl;
    return 1;
  }
  LenK   hC     { atof(argv[1]) };
  Mass   mL     { atof(argv[2]) };
  double alpha1 { atof(argv[3]) };

  ForceK thrust2Vac = 2.0 * 63'700.0_kg * g0K;
  ForceK thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  try
  {
    Ascent asc(mL, alpha1, thrust2Vac, thrust1Vac, &cout, 2);
    asc.FindFeasibleTrajectory(hC);
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
