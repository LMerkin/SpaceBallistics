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
#include <nlopt.hpp>

using namespace  SpaceBallistics;
namespace EAM  = EarthAtmosphereModel;
namespace LVAD = LVAeroDyn;

namespace
{
  //=========================================================================//
  // The "Ascent2" Class:                                                    //
  //=========================================================================//
  class Ascent2
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
    constexpr static decltype(Sqr(D)) CroS   = Pi<double> * Sqr(D) / 4.0;

    // ODE Integration Params: 1 msec step; it may only be reduced, never
    // increased beyond the original value:
    constexpr static Time   ODEInitStep      = 0.001_sec;
    constexpr static Time   ODEMaxStep       = ODEInitStep;
    constexpr static double ODERelPrec       = 1e-6;

    // Singular point detection criteria: NB: "Vhor" approaches 0 much faster
    // than "Vr":
    constexpr static VelK   SingVr           = VelK(1e-2); // 10   m/sec
    constexpr static VelK   SingVhor         = VelK(1e-4); //  0.1 m/sec

    // Angle-of-Attack (AoA) Limits for the 1st and the 2nd stage:
    constexpr static Angle  MaxAoA2          = To_Angle(10.0_deg);
    constexpr static Angle  MaxAoA1          = To_Angle( 1.0_deg);

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    // State Vector (for some time t <= 0, where t=t0=0 corresponds to the Orbi-
    // tal Insertion):
    // Computations are performed in TopoCentric (Start) Planar Polar CoOrds, so
    // the first 3 components are r, rDot, omega = phiDot. The 4th  component is
    // the Total Spent Propellant Mass between the curr "t" and "t0" (which is a
    // continuous and DECREASING function of "t",  as opposed to the Total Mass
    // which is discontinuous when the Stage1 or Fairing are jettisones), and
    // the 5th component is the polar angle "phi" (integrated omega):
    //
    using StateV  = std::tuple<LenK, VelK, AngVel, Mass, Angle>;
    //                         r   rDot=Vr omega   spent phi

    // The Time Derivative of the "StateV". The last component is the BurnRate:
    using DStateV = std::tuple<VelK, AccK, AngAcc, MassRate, AngVel>;

    // Flight Mode:
    enum class FlightMode: int
    {
      UNDEFINED = 0,
      Burn1     = 1,  // Stage1 Burn,
      Gap       = 2,  // Ballistic Gap
      Burn2     = 3   // Stage2 Burn
    };

    static char const* ToString(FlightMode a_mode)
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
      Angle const m_phi;
      Time  const m_t;

      // Non-Default Ctor:
      NearSingularityExn(LenK  a_r,   VelK a_vr, Mass a_spent_prop_mass,
                         Angle a_phi, Time a_t)
      : m_r             (a_r),
        m_Vr            (a_vr),
        m_spentPropMass (a_spent_prop_mass),
        m_phi           (a_phi),
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
    MassRate const        m_burnRateI2; // BurnRate @ Stage2 Ignition Time

    // Stage1 Params:
    Mass     const        m_fullMass1;
    Mass     const        m_emptyMass1;
    Mass     const        m_propMass1;
    ForceK   const        m_thrustVac1;
    MassRate const        m_burnRateI1; // BurnRate @ Stage1 Ignition Time

    // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
    // ignition):
    Time                  m_TGap;

    // Flight Control Program Parameterisation:
    //
    // BurnRate for Stages 1 and 2:
    // It is a non-increasing quadratic function of time, which coeffs depend-
    // ing (somehow) on 2 dimension-less params: "aHat" and "muHat":
    // BurnRate(tau) = BurnRateI + bMu * tau + aMu * tau^2, where "tau" is the
    // time since ignition of the corresp Stage:
    //
    using MassT2     =    decltype(MassRate(1.0) / 1.0_sec);
    using MassT3     =    decltype(MassT2  (1.0) / 1.0_sec);

    Time                  m_T2;         // Actual Stage2 BurnTime
    MassT3                m_aMu2;
    MassT2                m_bMu2;

    Time                  m_T1;         // Actual Stage1 BurnTime
    MassT3                m_aMu1;
    MassT2                m_bMu1;

    // AoA for Stage2: AoA(t) = t * (a * t + b), t <= 0:
    AngAcc                m_aAoA2;
    AngVel                m_bAoA2;

    // AoA for Stage1: AoA(nt)  = nt *  (a * nt  + b),
    // where nt = -tau =  tIgn1 - t <= 0;
    // thus,           AOA(tau) = tau * (a * tau - b):
    AngAcc                m_aAoA1;
    AngVel                m_bAoA1;

    // Transient Data (during flight path integration): St2 Cut-Off Time is 0:
    FlightMode            m_mode;
    Time                  m_ignTime2;       // Burn2 start: St2 Ignition Time
    Time                  m_fairingSepTime; // Fairing Separation Time
    Time                  m_cutOffTime1;    // Gap   start: St1 Cut-Off  Time
    Time                  m_ignTime1;       // Burn1 start: St1 Ignition Time

    // For output:
    std::ostream* const   m_os;
    int           const   m_logLevel;

  public:
    //=======================================================================//
    // Non-Default Ctor:                                                     //
    //=======================================================================//
    Ascent2
    (
      Mass          a_payload_mass,
      double        a_alpha1,      // FullMass1 / FullMass2
      ForceK        a_thrust2_vac,
      ForceK        a_thrust1_vac,
      std::ostream* a_os,
      int           a_log_level
    )
    : m_payLoadMass   (a_payload_mass),

      // Stage2 Params:
      m_fullMass2     ((StartMass - FairingMass - m_payLoadMass) /
                       (1.0 + a_alpha1)),
      m_emptyMass2    (m_fullMass2   * (1.0 - K2)),
      m_propMass2     (m_fullMass2   * K2),
      m_thrustVac2    (a_thrust2_vac),
      m_burnRateI2    (m_thrustVac2  / (IspVac2  * g0K)),

      // Stage1 Params:
      m_fullMass1     ((StartMass - FairingMass - m_payLoadMass) /
                       (1.0 + a_alpha1) * a_alpha1),
      m_emptyMass1    (m_fullMass1  * (1.0 - K1)),
      m_propMass1     (m_fullMass1  * K1),
      m_thrustVac1    (a_thrust1_vac),
      m_burnRateI1    (m_thrustVac1  / (IspVac1  * g0K)),
      m_TGap          (),     // 0 by default

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

      // Transient Data:
      m_mode          (FlightMode::UNDEFINED),
      m_ignTime2      (Time(NAN)),
      m_fairingSepTime(Time(NAN)),
      m_cutOffTime1   (Time(NAN)),
      m_ignTime1      (Time(NAN)),

      m_os            (a_os), // May be NULL
      m_logLevel      (a_log_level)
    {
      if ( IsNeg(m_payLoadMass) ||  IsNeg(m_TGap)      || !IsPos(m_fullMass2) ||
          !IsPos(m_thrustVac2)  || !IsPos(m_fullMass1) || !IsPos(m_thrustVac1))
        throw std::invalid_argument("Ascent2::Ctor: Invalid param(s)");

      // Therefore:
      assert(IsPos(m_emptyMass2) && IsPos(m_propMass2) && IsPos(m_burnRateI2) &&
             IsPos(m_emptyMass1) && IsPos(m_propMass1) && IsPos(m_burnRateI1));
    }

    //=======================================================================//
    // "RunRes":                                                             //
    //=======================================================================//
    // Possible results of "Run":
    //
    enum class RunRes: int
    {
      // The "desirable" result: Singular Point reached: V = 0 at some h >= 0:
      Singularity = 0,

      // Reached h = 0 at V > 0, still with unspent propellant (beyond the
      // minimum remnant):
      ZeroH       = 1,

      // Ran out of all available propellant (up to the minimum remnant), but
      // still h > 0 and V > 0:
      PropOut     = 2,

      // If any exception occurred:
      Error       = 3
    };

    static char const* ToString(RunRes a_rr)
    {
      switch (a_rr)
      {
        case RunRes::Singularity: return "Singularity";
        case RunRes::ZeroH      : return "ZeroH";
        case RunRes::PropOut    : return "PropOut";
        default                 : return "Error";
      }
    }

    //=======================================================================//
    // "Run": Integrate the Ascent Trajectory:                               //
    //=======================================================================//
    // Returns (RunRes, FinalH, FlightTime, ActStartMass):
    //
    std::tuple<RunRes, LenK, Time, Mass> Run
    (
      LenK    a_hc,

      // BurnRate ctl params for Stage2: The defaults corresp to const BurnRate:
      double  a_aHat2   = 0.0,     // Must be in [-1 .. 1]
      double  a_muHat2  = 1.0,     // Must be in [ 0 .. 1]

      // AoA ctl params for Stage2: The defaults corresp to AoA = 0:
      double  a_aAoA2   = 0.0,     // Must be in [ 0 .. 1]
      double  a_bAoA2   = 0.0,     // Must be in [ 0 .. 1]

      Time    a_t_gap   = 0.0_sec,

      // BurnRate ctl params for Stage1: The defaults corresp to const BurnRate:
      double  a_aHat1   = 0.0,     // Must be in [-1 .. 1]
      double  a_muHat1  = 1.0,     // Must be in [ 0 .. 1]

      // AoA ctl param for Stage1: The defaults corresp to AoA = 0:
      double  a_aAoA1   = 0.0,     // Must be in [ 0 .. 1]
      double  a_bAoA1   = 0.0      // Must be in [ 0 .. 1]
    )
    {
      //---------------------------------------------------------------------//
      // Check the Params:                                                   //
      //---------------------------------------------------------------------//
      if (a_hc < 100.0_km || IsNeg(a_t_gap))
        throw std::invalid_argument("Ascent2::Run: Invalid hC or TGap");

      // BurnRate ctls:
      if (std::fabs(a_aHat2) > 1.0 ||
          a_muHat2 <= 0.0          || a_muHat2 >= 1.0   ||
          a_aAoA2  <= 0.0          || a_aAoA2  >= 1.0   ||
          a_bAoA2  <= 0.0          || a_bAoA2  >= 1.0   ||
          std::fabs(a_aHat1) > 1.0 ||
          a_muHat1 <= 0.0          || a_muHat2 >= 1.0   ||
          a_aAoA1  <= 0.0          || a_aAoA1  >= 1.0   ||
          a_bAoA1  <= 0.0          || a_bAoA1  >= 1.0)
        throw std::invalid_argument
              ("Ascent2::Run: Invalid BurnRate or AoA param(s)");

      //---------------------------------------------------------------------//
      // BurnRate Coeffs:                                                    //
      //---------------------------------------------------------------------//
      // Normalise "aHat"s to [-1/3 .. +1/3]:
      a_aHat2 /= 3.0;
      a_aHat1 /= 3.0;

      // Compute the actual limits for "muHat", and the actual "muHat" within
      // those limits:
      double muHatLo2 = (1.0 - a_aHat2) / 2.0;
      double muHatUp2 = (a_aHat2 < 0.0) ? 1.0 + a_aHat2 : 1.0 - 2.0 * a_aHat2;
      assert(muHatLo2 <= muHatUp2);
      a_muHat2       *= (muHatUp2 - muHatLo2);
      // XXX: In a very degenerate case, we may get a_muHat2==0:
      assert(0.0 <= a_muHat2 && a_muHat2 <= 1.0);
      m_T2            = m_propMass2 / (m_burnRateI2 * a_muHat2);
      m_aMu2          = 3.0 * m_burnRateI2 / Sqr(m_T2) * a_aHat2;
      m_bMu2          = 2.0 * (m_propMass2 / Sqr(m_T2) - m_burnRateI2 / m_T2 -
                               m_aMu2 * m_T2 / 3.0);

      double muHatLo1 = (1.0 - a_aHat1) / 2.0;
      double muHatUp1 = (a_aHat1 < 0.0) ? 1.0 + a_aHat1 : 1.0 - 2.0 * a_aHat1;
      assert(muHatLo1 <= muHatUp1);
      a_muHat1 *= (muHatUp1 - muHatLo1);
      // XXX: In a very degenerate case, we may get a_muHat1==0:
      assert(0.0 <= a_muHat1 && a_muHat1 <= 1.0);
      m_T1            = m_propMass1 / (m_burnRateI1 * a_muHat1);
      m_aMu1          = 3.0 * m_burnRateI1 / Sqr(m_T1) * a_aHat1;
      m_bMu1          = 2.0 * (m_propMass1 / Sqr(m_T1) - m_burnRateI1 / m_T1 -
                               m_aMu1 * m_T1 / 3.0);

      //---------------------------------------------------------------------//
      // AoA Coeffs:                                                         //
      //---------------------------------------------------------------------//
      m_bAoA2        = - 4.0 * MaxAoA2 / m_T2 * a_bAoA2;
      AngAcc aAoALo2 = m_bAoA2 / m_T2;
      AngAcc aAoAUp2 =
        (a_bAoA2 < 0.5)
        ? (MaxAoA2 / m_T2 + m_bAoA2) / m_T2
        : - Sqr(m_bAoA2) / (4.0 * MaxAoA2);
      m_aAoA2 = aAoALo2  * (1.0 - a_aAoA2) + a_aAoA2 * aAoAUp2;

      m_bAoA1        = - 4.0 * MaxAoA1 / m_T1 * a_bAoA1;
      AngAcc aAoALo1 = m_bAoA1 / m_T1;
      AngAcc aAoAUp1 =
        (a_bAoA1 < 0.5)
        ? (MaxAoA1 / m_T1 + m_bAoA1) / m_T1
        : - Sqr(m_bAoA1) / (4.0 * MaxAoA1);
      m_aAoA1 = aAoALo1  * (1.0 - a_aAoA1) + a_aAoA1 * aAoAUp1;

      //---------------------------------------------------------------------//
      // For Testing Only:                                                   //
      //---------------------------------------------------------------------//
      if (m_os != nullptr && m_logLevel >= 2)
      {
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
      //---------------------------------------------------------------------//
      // Run the integration BACKWARDS from the orbital insertion point      //
      //---------------------------------------------------------------------//
      // (@ t=0):
      // Angular velocity at the Circular Orbit, XXX: with correction  for the
      // Earth Rotation Velocity (assuming the Target Orbit Inclination is the
      // same as the Launch Site Latitude):
      LenK   aC     = R + a_hc;
      VelK   v0     = SqRt(K / aC) - ERV;
      AngVel omega0 = 1.0_rad * v0 / aC;

      // The initial Radial Velocity is 0, and the final delta of Spent Mass
      // is also 0, and the initial Polar Angle is 0:
      StateV s0 = std::make_tuple(aC, VelK(0.0), omega0, 0.0_kg, 0.0_rad);

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
        //-------------------------------------------------------------------//
        // Actually Run the Integrator!                                      //
        //-------------------------------------------------------------------//
        tEnd =
          RKF5(&s0, t0, tMin, rhs,
               -ODEInitStep, -ODEMaxStep, ODERelPrec, &cb, m_os);
        assert(tEnd >=  tMin);
      }
      catch (NearSingularityExn const& ns)
      {
        //-------------------------------------------------------------------//
        // We have reached a vicinity of the singular point:                 //
        //-------------------------------------------------------------------//
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

          return std::make_tuple(RunRes::Singularity, singH, singT, singM);
        }
        else
          // Although we got a "NearSingularityExn", we could not determine the
          // precise location of the singular point:
          return std::make_tuple
                 (RunRes::Error, LenK(NAN), Time(NAN), Mass(NAN));
      }
      catch (std::exception const& exn)
      {
        //-------------------------------------------------------------------//
        // Any other "standard" (or other) exceptions: Ascent unsuccessful:  //
        //-------------------------------------------------------------------//
        if (m_os != nullptr && m_logLevel >= 1)
          *m_os  << "# Ascent2::Run: Exception: " << exn.what() << std::endl;
        return std::make_tuple(RunRes::Error, LenK(NAN), Time(NAN), Mass(NAN));
      }
      catch (...)
      {
        // Any other exception:
        if (m_os != nullptr && m_logLevel >= 1)
          *m_os  << "# Ascent2::Run: UnKnown Exception" << std::endl;
        return std::make_tuple(RunRes::Error, LenK(NAN), Time(NAN), Mass(NAN));
      }

      //---------------------------------------------------------------------//
      // Integration has run to completion, but not to the singular point:   //
      //---------------------------------------------------------------------//
      // Check the final state and mode:
      //
      LenK hEnd = std::get<0>(s0) - R;
      Mass mEnd = LVMass(s0, tEnd);

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
        return std::make_tuple(RunRes::PropOut, hEnd, tEnd, mEnd);
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
        return std::make_tuple(RunRes::ZeroH, 0.0_km, tEnd, mEnd);
      }
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

      // The result: NB: the last derivative is <= 0, since the corresp
      // "StateV" components is int_t^0 burnRate(t') dt', ie "t" is the
      // LOWER integration limit:
      return std::make_tuple(Vr, r2Dot, omegaDot, -burnRate, omega);
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
      Angle  phi           = std::get<4>(*a_s);

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
        throw NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

      //---------------------------------------------------------------------//
      // Generic Case:                                                       //
      //---------------------------------------------------------------------//
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
      //---------------------------------------------------------------------//
      // Switching Gap -> Burn1:                                             //
      //---------------------------------------------------------------------//
      // Occurs merely by time delay (if the Gap is 0, this will happen immedi-
      // ately after the above change Burn2 -> Gap):
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

    //=======================================================================//
    // Atmospheric Conditions and Aerodynamic Drag Force:                    //
    //=======================================================================//
    static std::tuple<EAM::AtmConds, ForceK, ForceK>
    AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA)
    {
      // Curr altitude (possible slightly negative vals are rounded to 0):
      auto atm  = EAM::GetAtmConds(a_r - R);

      // The Mach Number:
      double M  = double(To_Len_m(a_v) / std::get<3>(atm));
      assert(M >= 0.0);

      // The Aerodynamic Force Main Term:
      ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * CroS;

      // The Drag and Lift Coeffs:
      double cD = LVAD::cD(M, a_AoA);
      double cL = LVAD::cL(M, a_AoA);

      return std::make_tuple(atm, cD * F, cL * F);
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
        // Then the LV is falling back to the pad, and the singular point is
        // not reachable; this is because we have arrived at the mass "m1"
        // which is too large:
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

  public:
    //=======================================================================//
    // Propellant Burn Rate: May be variable:                                //
    //=======================================================================//
    MassRate PropBurnRate(Time a_t) const
    {
      switch (m_mode)
      {
      case FlightMode::Burn1:
      {
        // Arg: "t" is the time from Stage1 cut-off (so -T1 <= t <= 0):
        Time t = std::max(std::min(a_t - m_cutOffTime1, 0.0_sec), -m_T1);
        return   std::max(m_burnRateI1 + (m_bMu1 + m_aMu1 * t) * t, 
                          MassRate(0.0));
      }
      case FlightMode::Gap:
        return MassRate(0.0);

      case FlightMode::Burn2:
      {
        // Arg: "t" is time from Stage2 cut-off (same as orbital insertion
        // time), so it is just the main time "a_t": -T2 <= t <= 0; enforce
        // those constraints:
        assert(!IsPos(a_t));
        Time t = std::max(a_t, -m_T2);
        return   std::max(m_burnRateI2 + (m_bMu2 + m_aMu2 * t) * t,
                          MassRate(0.0));
      }
      default:
        return MassRate(0.0);
      }
    }

    //=======================================================================//
    // Angle-of-Attack:                                                      //
    //=======================================================================//
    Angle AoA(Time a_t) const
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
        return std::max(tau * (m_aAoA1 * tau - m_bAoA1), 0.0_rad);
      }
      case FlightMode::Gap:
        return 0.0_rad;

      case FlightMode::Burn2:
      {
        // Arg: "t" is time from Stage2 cut-off (same as orbital insertion
        // time), so it is just the main time "a_t": -T2 <= t <= 0; enforce
        // those constraints:
        assert(!IsPos(a_t));
        Time t = std::max(a_t, -m_T2);
        return   std::max(t *  (m_aAoA2 * t + m_bAoA2), 0.0_rad);
      }
      default:
        return 0.0_rad;
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
    // "FindFeasibleTGap"                                                    //
    //=======================================================================//
    // "TGap" is a control param  which must ensure a correct ascent to the
    // "a_hc" circular orbit. This method tries to find a suitable value of
    // "TGap" such that  the singular point of the ascent trajectory exists
    // and is @ h=0. Returns "TGap" if found, or nullopt otherwise:
    //
    std::optional<Time> FindFeasibleTGap(LenK a_hc)
    {
      constexpr Time MaxTGap = 300.0_sec;
      Time           tgStep  = 10.0_sec;
      Time           tg0     = - tgStep;  // As no valid "tg0" is found yet
      LenK           h0(NAN);             // Not known yet
      LenK           MaxH0   = 0.5_km;    // Max acceptable h0

      // Run iterations for "tg" multiples times and to find the "tg" corresp
      // to h ~= 0:
      for (Time tg = tg0 + tgStep; tg <= MaxTGap; tg += tgStep)
      {
        if (IsNeg(tg))
          // This may happen if no Singular Points were found at all:
          return std::nullopt;

        //-------------------------------------------------------------------//
        // Run the integrator and hope to reach the singular point:          //
        //-------------------------------------------------------------------//
        auto res =
          Run
          (
            // Circular Orbit Altitude:
            a_hc,

            // Stage2 BurnRate Ctls:
            -0.5, // aHat2
            0.5,  // muHat2

            // Stage2 AoA Ctls:
            0.2,  // aAoA2
            0.1,  // bAoA2

            // Ballistic Gap:
            tg,

            // Stage1 BurnRate Ctls:
            0.3,  // aHat1
            0.7,  // muHat1

            // Stage1 AoA Ctls:
            0.2,  // aAoA1
            0.9   // bAoA1
          );
        //-------------------------------------------------------------------//
        // Analyse the results:                                              //
        //-------------------------------------------------------------------//
        RunRes rr  = std::get<0>(res);
        LenK   h   = std::get<1>(res);
        Time   t   = std::get<2>(res);
        Mass   m   = std::get<3>(res);

        if (m_os != nullptr)
          *m_os << "# t=" << t.Magnitude()     << " sec, h="
                << h.Magnitude()  << " km, m=" << m.Magnitude() << " kg: "
                << ToString(rr)   << " w/ TGap="
                << tg.Magnitude() << std::endl << std::endl;

        switch (rr)
        {
        case RunRes::Singularity:
          // If there was already a valid "h0", a new valid one ("h") should
          // always be smaller than "h0":
          if (IsFinite(h0) && h > h0 && m_os != nullptr)
            *m_os << "Ascent2::FindFeasibleTGap: Non-Monotonic Singular Point: "
                     "tg0="      << tg0.Magnitude()
                  << " sec, h0=" << h0.Magnitude()
                  << " km; tg="  << tg.Magnitude()
                  << " sec; h="  << h.Magnitude()  << std::endl;

          // Yet in any case, update "tg0" which is the maximum "tg" corresp to
          // the Singularity (therefore with the minimum "h") found yet:
          tg0 = tg;
          h0  = h;
          // Continue iterations:
          break;

        case RunRes::ZeroH:
          // There will be no Singularities for larger "tg"s; so either reduce
          // the step or accept the memoised "tg0" and "h0":
          if (tgStep >= 1.0_sec)
          {
            tgStep /= 10.0;
            tg      = tg0;   // Will continue with tg0 + tgStep_new
            continue;
          }
          else
          {
            // The "tgStep" is already small enough, so accept the last valid
            // "tg0". Check that the corresp "h0" is sufficiently small:
            assert(!IsNeg(h0));
            if (h0 > MaxH0 && m_os != nullptr)
              *m_os << "Ascent2::FindFeasibleTGap: TGap="    << tg0.Magnitude()
                    << " sec: Imprecise Singular Point: h=" << h0.Magnitude()
                    << " km" << std::endl;
            return tg0;
          }

        default:
          // "PropOut" or "Error":
          // There will be no Singularities anymore, and it is unlikely that we
          // would get h ~= 0 for any "tg". However, if we already got a more
          // or less acceptable one, return it:
          if (h0 <= MaxH0)
            return tg0;
          else
            return std::nullopt;
        }
      }
      // If we got here: We have iterated to "MaxTGap" without success:
      return std::nullopt;
    }
  };
  // End of "Ascent2" Class
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

  ForceK thrust2Vac = 1.0 * 63'700.0_kg * g0K;
  ForceK thrust1Vac = 9.0 * 59'500.0_kg * g0K;
  try
  {
    Ascent2   asc(mL, alpha1, thrust2Vac, thrust1Vac, &cout, 2);

    auto tg = asc.FindFeasibleTGap(hC);

    if (bool(tg))
      std::cout << tg.value().Magnitude() << " sec" << std::endl;
    else
      std::cout << "None" << std::endl;
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
	return  0;
}
