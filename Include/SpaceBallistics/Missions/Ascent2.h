// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/Missions/Ascent2.h":                 //
//                   Ascent-to-Orbit for a "Model" 2-Stage LV                //
//===========================================================================//
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.hpp"
#include <boost/functional/hash.hpp>
#include <ostream>
#include <optional>
#include <utility>
#include <unordered_map>

namespace SpaceBallistics
{
  namespace EAM = EarthAtmosphereModel;

  //=========================================================================//
  // The "Ascent2" Class:                                                    //
  //=========================================================================//
  class Ascent2
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // TODO: Some of those Consts should be made dynamically-configurable par-
    // ams:
    constexpr static GMK    K                = BodyData<Body::Earth>::K;
    constexpr static LenK   R                = 6371.0_km;   // Equi-Volume

    // LV Params:
    // Stage1 (engine performance @ SL and in Vac); its full mass is determined
    // dynamically, but we need some initial estimates for reference:
    constexpr static double K1               = 0.9;
//  constexpr static Time   IspSL1           = 326.0_sec;
    constexpr static Time   IspSL1           = 320.0_sec;
//  constexpr static Time   IspVac1          = 353.0_sec;
    constexpr static Time   IspVac1          = 346.6_sec;
    constexpr static double PropRem1         = 0.01;

    // Stage2 (only in Vac): Most params are dynamic. Here "PropRem2" is a frac-
    // tion of the total Stage2 propellant load, not a fixed value:
    constexpr static double K2               = 0.93;
//  constexpr static Time   IspVac2          = 381.0_sec;
    constexpr static Time   IspVac2          = 377.0_sec;
    constexpr static double PropRem2         = 0.01;

    // Total Start Mass (Constant):
    constexpr static Mass   MaxStartMass     = 360'000.0_kg;

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
    constexpr static Angle  MaxAoA2          = To_Angle(20.0_deg);
    constexpr static Angle  MaxAoA1          = To_Angle( 2.0_deg);

    // The Deepest Throttling Level (ie the Min BurnRate relative to the Max
    // one) for all engines:
    constexpr static double DT               = 0.5;

    // The Number of Ctl Params in the optimisation algorithm for construction
    // of the ascent-to-orbit trajectory:
    // [bHat2, muHat2, aAoA2, bAoA2, TGap, bHat1, muHat1, aAoA1, bAoA1], so:
    constexpr static int    NP               = 9;

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // "StateV":                                                             //
    //-----------------------------------------------------------------------//
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

    //-----------------------------------------------------------------------//
    // "DStateV":                                                            //
    //-----------------------------------------------------------------------//
    // The Time Derivative of the "StateV". The "MassRate" components is the
    // negated (<= 0) BurnRate:
    //
    using DStateV = std::tuple<VelK, AccK, AngAcc, MassRate, AngVel>;

    //-----------------------------------------------------------------------//
    // "FlightMode":                                                         //
    //-----------------------------------------------------------------------//
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
    //-----------------------------------------------------------------------//
    // "NearSingularityException" Class:                                     //
    //-----------------------------------------------------------------------//
    // Exception thrown when the flight path is approaching the singular point:
    //
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

    //-----------------------------------------------------------------------//
    // "RunRC":                                                              //
    //-----------------------------------------------------------------------//
    // Possible results of "Run":
    //
    enum class RunRC: int
    {
      // The "desirable" result: Singular Point reached: V = 0 at some h >= 0:
      Singularity = 0,

      // Reached h = 0 at V > 0, still with unspent propellant (beyond the
      // minimum remnant):
      ZeroH       = 1,

      // Ran out of all available propellant (up to the minimum remnant), but
      // still h > 0 and V > 0:
      FlameOut    = 2,

      // If any exception occurred:
      Error       = 3
    };

    static char const* ToString(RunRC a_rc)
    {
      switch (a_rc)
      {
        case RunRC::Singularity: return "Singularity";
        case RunRC::ZeroH      : return "ZeroH";
        case RunRC::FlameOut   : return "FlameOut";
        default                : return "Error";
      }
    }

    //-----------------------------------------------------------------------//
    // "RunRes" Struct:                                                      //
    //-----------------------------------------------------------------------//
    struct RunRes
    {
      RunRC m_rc;   // Return Code
      Time  m_T;    // Final (actually start) Time, < 0
      LenK  m_hT;   // Final (actually start) Altitude
      VelK  m_VT;   // Final (sctually start) Velocity
      Mass  m_mT;   // Final (actually start) Mass
    };

  public:
    //-----------------------------------------------------------------------//
    // "AscCtlsL" Struct:                                                    //
    //-----------------------------------------------------------------------//
    // Dimension-Less Params Controlling the Ascent Trajectory:
    //
    struct AscCtlsL
    {
      //---------------------------------------------------------------------//
      // Data Flds:                                                          //
      //---------------------------------------------------------------------//
      // BurnRate Ctl Params for Stage2: The defaults corresp to const BurnRate:
      double  m_bHat2;            // Must be in [0 .. 1], default is 0
      double  m_muHat2;           // Must be in [0 .. 1], default is 1

      // AoA Ctl Params for Stage2: The defaults corresp to AoA = 0:
      double  m_aAoA2;            // Must be in [0 .. 1], default is 0
      double  m_bAoA2;            // Must be in [0 .. 1], default is 0

      // Ballistic Gap: XXX: Using "double" rather than "Time" here, for STL
      // compatibility:
      double  m_TGap;             // Default is 0

      // BurnRate ctl params for Stage1: The defaults corresp to const BurnRate:
      double  m_bHat1;            // Must be in [0 .. 1], default is 0
      double  m_muHat1;           // Must be in [0 .. 1], default is 1

      // AoA ctl param for Stage1: The defaults corresp to AoA = 0:
      double  m_aAoA1;            // Must be in [0 .. 1], default is 0
      double  m_bAoA1;            // Must be in [0 .. 1], default is 0

      //---------------------------------------------------------------------//
      // To use "AscCtlsL" as an "unordered_map" key, we need Hash and (==): //
      //---------------------------------------------------------------------//
      bool operator==(AscCtlsL const&) const = default;

      // For Hash, see the end of this header...
    };

    //-----------------------------------------------------------------------//
    // "AscCtlsD": Dimensioned Params Controlling the Ascent Trajectory:     //
    //-----------------------------------------------------------------------//
    using MassT2 = decltype(MassRate(1.0) / 1.0_sec);
    using MassT3 = decltype(MassT2  (1.0) / 1.0_sec);

    struct AscCtlsD
    {
      // Similar to the corresp flds in "Ascent2" class itself:
      // BurnRate for Stage2:
      Time                  m_T2;         // Actual Stage2 BurnTime
      MassT3                m_aMu2;
      MassT2                m_bMu2;
      // AoA for Stage2: AoA(t) = t * (a * t + b), t <= 0:
      AngAcc                m_aAoA2;
      AngVel                m_bAoA2;

      // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
      // ignition):
      Time                  m_TGap;

      // BurnRate for Stage1:
      Time                  m_T1;         // Actual Stage1 BurnTime
      MassT3                m_aMu1;
      MassT2                m_bMu1;

      // AoA for Stage1: AoA(nt)  = nt *  (a * nt  + b),
      // where nt = -tau =  tIgn1 - t <= 0;
      // thus,           AOA(tau) = tau * (a * tau - b):
      AngAcc                m_aAoA1;
      AngVel                m_bAoA1;

      // In addition, we provide the LV StartMass (<= MaxStartMass):
      Mass                  m_startMass;
    };

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // LV Params                                                             //
    //-----------------------------------------------------------------------//
    // XXX: In the future, these params may become "constexpr"s, or other way
    // round,  the params which are now "constexpr"s may be moved to here for
    // the maximum flexibility...
    //
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

    //-----------------------------------------------------------------------//
    // Mission Params:                                                       //
    //-----------------------------------------------------------------------//
    Mass     const        m_payLoadMass;
    LenK                  m_Rins;       // Radius-Vector @ Orbital Insertion
    VelK                  m_Vins;       // LV Velocity   @ Orbital Insertion

    //-----------------------------------------------------------------------//
    // Flight Control Program Parameterisation:                              //
    //-----------------------------------------------------------------------//
    // BurnRate for Stage2:
    // It is a non-increasing quadratic function of time:
    // BurnRate(tau) = BurnRateI + bMu * tau + aMu * tau^2, where "tau" is the
    // time since Stage2 ignition (0 <= tau <= T2):
    //
    Time                  m_T2;         // Actual Stage2 BurnTime
    MassT3                m_aMu2;
    MassT2                m_bMu2;

    // AoA for Stage2: AoA(t) = t * (a * t - b),
    // where t <= 0 is our standard integration backward-running time    (t=0
    // corresponds  to the orbital insertion instant, and AoA=0 at that point);
    // NB: the (-b) coeff !!!
    AngAcc                m_aAoA2;
    AngVel                m_bAoA2;

    // Ballistic Gap (a passive interval between Stage1 cut-off and Stage2
    // ignition):
    Time                  m_TGap;

    // BurnRate for Stage1: Similar to that of Stage1, where "tau" is the time
    // since Stage1 ignition (0 <= tau <= T1):
    Time                  m_T1;         // Actual Stage1 BurnTime
    MassT3                m_aMu1;
    MassT2                m_bMu1;

    // AoA for Stage1: for similarity with Stage2,
    // AoA(nt)  = nt * (a * nt  - b),
    // where nt = -tau  = tIgn1 - t <= 0,  so AoA=0 @ nt=0 (Stage1 ignition);
    // thus,   AOA(tau) = tau * (a * tau + b),
    // where "tau" is again the time since Stage1 ignition (0 <= tau <= T1):
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

    // Static Cache for "RunRes"es (for use with NLOpt):
    static std::unordered_map<AscCtlsL, RunRes> s_Cache;

  public:
    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    Ascent2
    (
      // LV Params (XXX: they should either be all "constexpr"s or all configu-
      // rable; but for now, some params are "fixed" and some are not):
      double          a_alpha1,         // FullMass1 / FullMass2
      ForceK          a_thrust2_vac,
      ForceK          a_thrust1_vac,

      // Mission Params:
      Mass            a_payload_mass,
      LenK            a_h_perigee,
      LenK            a_h_apogee,
      Angle_deg       a_incl,
      Angle_deg       a_launch_lat,

      // Logging Params:
      std::ostream*   a_os,
      int             a_log_level
    );

    //-----------------------------------------------------------------------//
    // Copy Ctor: Required for NLOpt-based optimisations:                    //
    //-----------------------------------------------------------------------//
    Ascent2(Ascent2 const&) = default;

    //-----------------------------------------------------------------------//
    // "Run": Integrate the Ascent Trajectory:                               //
    //-----------------------------------------------------------------------//
    RunRes Run();

  private:
    //-----------------------------------------------------------------------//
    // Internal Helpers:                                                     //
    //-----------------------------------------------------------------------//
    // "SetCtlsParams":
    // Setting the Flight Ctl Params. Used in the Ctor and in the optimisation
    // loop:
    void SetCtlParams(AscCtlsL const& a_ctls);

    void SetCtlParams
    (
      double a_bHat2, double a_muHat2,  // Stage2 BurnRate
      double a_aAoA2, double a_bAoA2,   // Stage2 AoA
      double a_TGap,                    // Ballistic Gap in sec
      double a_bHat1, double a_muHat1,  // Stage1 BurnRate
      double a_aAoA1, double a_bAoA1    // Stage1 AoA
    );

    // "ODERHS":
    // For Ascent Trajectory Integration:
    //
    DStateV ODERHS(StateV const& a_s, Time a_t);

    // "ODECB":
    // CallBack (invoked after the completion of each integration step):
    // Here FlightMode switching occurs, so this method is non-"const":
    //
    bool ODECB(StateV* a_s, Time a_t);

    // "AeroDynForces":
    // Atmospheric Conditions and Aerodynamic Drag Force:
    //
    static std::tuple<EAM::AtmConds, ForceK, ForceK>
    AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA);

    // "LocateSingularPoint":
    // The final stage of integration, where we assume omega=0 (a purely vert-
    // ical motion)   and integrate the simplified ODEs ANALYTICALLY to avoid
    // numerical instabilities in the viciniy of the singular point.
    // Returns (singH, singT)  if the singular point has been found, otherwise
    // "nullopt":
    //
    std::optional<std::pair<StateV, Time>>
    LocateSingularPoint(NearSingularityExn const& a_nse);

    // "PropBurnRate": (may be variable over time, >= 0):
    //
    MassRate PropBurnRate(Time a_t) const;

    // "AoA": Angle-of-Attack (variable over time):
    //
    Angle AoA(Time a_t) const;

    // "Thrust": Depends on the Mode, BurnRate and the Counter-Pressure:
    //
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const;

    // "LVMass": Current Mass (LV + PayLoad):
    //
    Mass LVMass(StateV const& a_s, Time a_t) const;

    // "OutputCtls": For Testing Only:
    //
    void OutputCtls() const;

    //-----------------------------------------------------------------------//
    // NLOpt Support:                                                        //
    //-----------------------------------------------------------------------//
    static std::optional<Ascent2::RunRes> GetRunRes
           (double const a_xs[NP], void* a_env);

    static double EvalNLOptObjective
    (
      unsigned     a_n,
      double const a_xs[NP],
      double*      a_grad,
      void*        a_env
    );

    static void EvalNLOptConstraints
    (
      unsigned     a_m,
      double       a_constrs[],
      unsigned     a_n,
      double const a_xs[NP],
      double*      a_grad,
      void*        a_env
    );

  public:
    //-----------------------------------------------------------------------//
    // "FindOptimalAscentCtls"                                               //
    //-----------------------------------------------------------------------//
    // Tries to find the Flight Control Params such that  the specified target
    // orbit is reached (with the orbital insertion occurring in the perigee,
    // if the orbit is elliptical) and that requires the minimum LV start mass
    // (whereas the payload mass is fixed).
    // Returns (OptAscCtlsD, MinStartMass):
    //
    static std::optional<AscCtlsD> FindOptimalAscentCtls
    (
      // LV Params (those which are considered to be non-"constexpr"):
      double          a_alpha1,      // FullMass1 / FullMass2
      ForceK          a_thrust2_vac,
      ForceK          a_thrust1_vac,

      // Mission Params:
      Mass            a_payload_mass,
      LenK            a_h_perigee,
      LenK            a_h_apogee,
      Angle_deg       a_incl,
      Angle_deg       a_launch_lat,

      // Logging Params:
      std::ostream*   a_os,
      int             a_log_level
    );
  };
}
// End namespace SpaceBallistics

namespace std
{
  //=========================================================================//
  // Specialize "std::hash" for "Ascent2::AscCtlsL":                         //
  //=========================================================================//
  // Using "boost::hash_combine":
  //
  template <>
  struct hash<SpaceBallistics::Ascent2::AscCtlsL>
  {
    size_t operator() (SpaceBallistics::Ascent2::AscCtlsL const& c) const
    {
      // Normalize -0.0 to 0.0 to avoid inconsistency with (==):
      auto norm0 = [](double a_x) -> double
                   { return (a_x  == 0.0) ? 0.0 : a_x; };

      size_t seed = 0;
      boost::hash_combine(seed, norm0(c.m_bHat2 ));
      boost::hash_combine(seed, norm0(c.m_muHat2));
      boost::hash_combine(seed, norm0(c.m_aAoA2 ));
      boost::hash_combine(seed, norm0(c.m_bAoA2 ));
      boost::hash_combine(seed, norm0(c.m_TGap  ));
      boost::hash_combine(seed, norm0(c.m_bHat1 ));
      boost::hash_combine(seed, norm0(c.m_muHat1));
      boost::hash_combine(seed, norm0(c.m_aAoA1 ));
      boost::hash_combine(seed, norm0(c.m_bAoA1 ));
      return seed;
    }
  };
}
// End namespace "std"
