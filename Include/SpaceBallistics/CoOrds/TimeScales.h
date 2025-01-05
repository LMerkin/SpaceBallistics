// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/TimeScales.h":                  //
//===========================================================================//
#pragma   once
#include "SpaceBallistics/Types.hpp"
#include <tuple>

namespace SpaceBallistics
{
  //=========================================================================//
  // Fwd Decls:                                                              //
  //=========================================================================//
  class TT;
  class TDB;

  //=========================================================================//
  // Consts:                                                                 //
  //=========================================================================//
  // The J2000.0 Epoch = 2000 Jan 1.5 (will be used in TT and TDB):
  //
  constexpr inline Time_day Epoch_J2000    = 2451545.0_day;
  constexpr inline Time_jyr Epoch_J2000_Yr = Time_jyr(2000.0 + 0.5 / 365.25);

  //=========================================================================//
  // "UTC" Struct:                                                           //
  //=========================================================================//
  // Terrestrial Non-Uniform Civil Time. For external representation of  time
  // only; no computations are performed directly with "UTC";  it needs to be
  // converted to "TT" first, and then possibly to "TDB".
  // XXX: The user can create any (even invalid) "UTC" struct as a mere aggre-
  // gate, but it will be validated when converted to "TT":
  //
  struct UTC
  {
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    int    const m_year;    //  Gregorian Calendar Year, assumed to be > 0
    int    const m_month;   //  1..12 (NOT using the UNIX convention, 0..11!)
    int    const m_day;     //  1..DaysInMonth(Year, Month)
    int    const m_hour;    //  0..23
    int    const m_min;     //  0..59
    double const m_sec;     // [0..60), or [0..61) if LeapSecond

    // NB: Default Ctor is implicitly deleted!

    //-----------------------------------------------------------------------//
    // Non-Default Ctors: From Date and Time:                                //
    //-----------------------------------------------------------------------//
    constexpr explicit UTC
    (
      std::tuple<int, int, int>    const& a_date,
      std::tuple<int, int, double> const& a_time = std::make_tuple(0, 0, 0.0)
    );

    constexpr explicit UTC
    (
      int a_year,     int a_month = 1, int    a_day = 1,
      int a_hour = 0, int a_min   = 0, double a_sec = 0.0
    );

    //-----------------------------------------------------------------------//
    // Comparison:                                                           //
    //-----------------------------------------------------------------------//
    constexpr bool operator== (UTC const& a_r) const = default;
    constexpr bool operator<  (UTC const& a_r) const;

    //-----------------------------------------------------------------------//
    // Leap Years and Leap Seconds:                                          //
    //-----------------------------------------------------------------------//
    // "IsLeapYear":
    // Is this year a Gregorian Leap Year?
    //
    constexpr static bool IsLeapYear(int a_year);

    // "DaysInMonth":
    constexpr static int DaysInMonth(int a_year, int a_month);

    // "HasLeapSecond":
    // Whether this "UTC" instant (considered up to a Minute) contains a Leap
    // Second:
    //
    constexpr bool HasLeapSecond() const;

    // "GetCumPrevLeapSecs":
    // Cumulative number of Leap Seconds BEFORE this instant. If it corresponds
    // to a Leap Second itself, the current Leap Second is NOT counted yet:
    //
    constexpr Time GetCumPrevLeapSecs() const;

    //-----------------------------------------------------------------------//
    // "IsValid": Over-All UTC Validation:                                   //
    //-----------------------------------------------------------------------//
    constexpr bool IsValid() const;
  };

  //=========================================================================//
  // "TT" Class:                                                             //
  //=========================================================================//
  // Terrestrial Time (Uniform, Relativistic), implemented via TAI (with an
  // offset):
  //
  class TT
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // The time is stored in "Modified Julian Seconds",  ie as a number of TAI
    // seconds since the (Modified) Epoch of J2000.0. XXX: This representation
    // has precision of 1 microsecond over the time span of +- 150 years only,
    // which is barely sufficient:
    Time      m_MJS;

  public:
    //-----------------------------------------------------------------------//
    // Default Ctor, Copy Ctor, Assignment, Equality:                        //
    //-----------------------------------------------------------------------//
    // NB: Default Ctor actually returns the Epoch_J2000 in TT:
    constexpr TT              ()                = default;
    constexpr TT              (TT const&)       = default;
    constexpr TT&  operator=  (TT const&)       = default;
    constexpr bool operator== (TT const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // From UTC:
    constexpr explicit  TT(UTC const& a_utc);

    // Directly constructing TT from JD_TT:
    constexpr explicit TT(Time_day a_jd_tt)
    : m_MJS    (To_Time_sec(a_jd_tt - Epoch_J2000))
    {}

    // Directly constructing from JYear_TT:
    constexpr explicit TT(Time_jyr a_jyr_tt)
    : m_MJS    (To_Time_sec( a_jyr_tt - Epoch_J2000_Yr))
    {}

    // Directly constructing from MJS_TT: USE WITH CARE!
    constexpr explicit TT(Time a_mjs_tt)
    : m_MJS    (a_mjs_tt)
    {}

    // From TDB (not "constexpr", as DE440T is required):
    explicit TT(TDB const& a_tdb);

    // Constructing an "UnDef":
    constexpr static TT UnDef()       { return TT(Time(NaN<double>)); }

    // Checking for an "UnDef":
    constexpr bool    IsUnDef() const { return IsNaN(m_MJS); }

    //-----------------------------------------------------------------------//
    // Extracting the Time value (since the Epoch) or JD_TT val:             //
    //-----------------------------------------------------------------------//
    // XXX: USE WITH CARE, as TimeScale info is then lost:
    //
    constexpr Time     GetTimeSinceEpoch() const
      { return m_MJS; }

    constexpr Time_day GetJDsSinceEpoch () const
      { return To_Time_day(m_MJS); }

    // Extracting the JD_TT:
    constexpr Time_day GetJD            () const
      { return Epoch_J2000    + To_Time_day(m_MJS); }

    // Extracting the JYr_TT:
    constexpr Time_jyr GetJYr           () const
      { return Epoch_J2000_Yr + To_Time_jyr(m_MJS); }

    //-----------------------------------------------------------------------//
    // Time Intervals (Durations):                                           //
    //-----------------------------------------------------------------------//
    constexpr Time operator-  (TT   a_right) const
      { return m_MJS - a_right.m_MJS;    }

    constexpr TT&  operator+= (Time a_right)
    {
      m_MJS += a_right;
      return  *this;
    }

    constexpr TT&  operator-= (Time a_right)
    {
      m_MJS -= a_right;
      return  *this;
    }

    constexpr TT   operator+  (Time a_right) const
    {
      TT res = *this;
      res   += a_right;
      return   res;
    }

    constexpr TT   operator-  (Time a_right) const
    {
      TT res = *this;
      res   -= a_right;
      return   res;
    }

    //-----------------------------------------------------------------------//
    // Comparisons:                                                          //
    //-----------------------------------------------------------------------//
    constexpr bool operator<  (TT a_right) const
      { return m_MJS <  a_right.m_MJS; }

    constexpr bool operator<= (TT a_right) const
      { return m_MJS <= a_right.m_MJS; }

    constexpr bool operator>  (TT a_right) const
      { return m_MJS >  a_right.m_MJS; }

    constexpr bool operator>= (TT a_right) const
      { return m_MJS >= a_right.m_MJS; }

    constexpr bool ApproxEquals
      (TT a_right, double a_tol = DefaultTol<double>) const
      { return m_MJS.ApproxEquals(a_right.m_MJS, a_tol); }

  public:
    //-----------------------------------------------------------------------//
    // TT -> UTC Conversion:                                                 //
    //-----------------------------------------------------------------------//
    constexpr operator UTC() const;
  };

  //=========================================================================//
  // "TDB" Class:                                                            //
  //=========================================================================//
  // Solar System BaryCentric Time (Uniform, Relativistic). Used as the argument
  // in JPL Ephemerides (DE440T):
  //
  class TDB
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // The time is stored in "Julian Seconds", ie as a number of TAI seconds
    // since the Epoch.  XXX: Again, this representation  has  precision  of
    // 1 microsecond over the time span of +- 150 years only, which is barely
    // sufficient:
    Time  m_MJS;

  public:
    //-----------------------------------------------------------------------//
    // Default Ctor, Copy Ctor, Assignment, Equality:                        //
    //-----------------------------------------------------------------------//
    // NB: Default Ctor returns the Epoch_J2000 but in TDB:
    constexpr TDB              ()                 = default;
    constexpr TDB              (TDB const&)       = default;
    constexpr TDB&  operator=  (TDB const&)       = default;
    constexpr bool  operator== (TDB const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // Constructing TDB from TT. The implementation is non-trivial, requires
    // the JPL Ephemerides (DE440T), so it is not a "constexpr":
    explicit TDB(TT const& a_tt);

    // Directly constructing TDB from JD_TDB:
    constexpr explicit TDB (Time_day a_jd_tdb)
    : m_MJS     (To_Time_sec(a_jd_tdb - Epoch_J2000))
    {}

    // Directly constructing TDB from MJS_TDB: USE WITH CARE!
    constexpr explicit TDB (Time a_mjs_tdb)
    : m_MJS     (a_mjs_tdb)
    {}

    // Constructing an  "UnDef":
    constexpr static TDB UnDef()       { return TDB(Time(NaN<double>)); }

    // Checking for an  "UnDef":
    constexpr bool     IsUnDef() const { return IsNaN(m_MJS); }

    //-----------------------------------------------------------------------//
    // Extracting the Time value (since the Epoch_J2000) or JD_TDB val:      //
    //-----------------------------------------------------------------------//
    // XXX: USE WITH CARE, as TimeScale info is then lost:
    //
    constexpr Time     GetTimeSinceEpoch() const
      { return m_MJS; }

    constexpr Time_day GetJDsSinceEpoch () const
      { return To_Time_day(m_MJS); }

    // Extracting the JD_TDB:
    constexpr  Time_day GetJD           () const
      { return Epoch_J2000 + To_Time_day(m_MJS); }

    // XXX: Unlike TT, for TDB we probably don't need "GetJYr"...

    //-----------------------------------------------------------------------//
    // Time Intervals (Durations):                                           //
    //-----------------------------------------------------------------------//
    constexpr Time operator-  (TDB  a_right) const
      { return m_MJS - a_right.m_MJS;    }
    
    constexpr TDB& operator+= (Time a_right)
    {
      m_MJS += a_right;
      return  *this;
    }

    constexpr TDB& operator-= (Time a_right)
    {
      m_MJS -= a_right;
      return  *this;
    }

    constexpr TDB  operator+  (Time a_right) const
    {
      TDB res = *this;
      res    += a_right;
      return    res;
    }

    constexpr TDB  operator-  (Time a_right) const
    {
      TDB res = *this;
      res    -= a_right;
      return    res;
    }

    //-----------------------------------------------------------------------//
    // Comparisons:                                                          //
    //-----------------------------------------------------------------------//
    constexpr bool operator<  (TDB a_right) const
      { return m_MJS <  a_right.m_MJS; }

    constexpr bool operator<= (TDB a_right) const
      { return m_MJS <= a_right.m_MJS; }

    constexpr bool operator>  (TDB a_right) const
      { return m_MJS >  a_right.m_MJS; }

    constexpr bool operator>= (TDB a_right) const
      { return m_MJS >= a_right.m_MJS; }

    constexpr bool  ApproxEquals
      (TDB a_right, double a_tol = DefaultTol<double>) const
      { return m_MJS.ApproxEquals(a_right.m_MJS, a_tol); }
  };

  //=========================================================================//
  // "FlightTime":                                                           //
  //=========================================================================//
  // In seconds since LaunchTime (given as TT). XXX: Unfortunately, Launch Time
  // cannot be provided statically; and using TT might not be a good idea for
  // launches from extraterrestrial bodies; but for now, we are primarily con-
  // cerned with launches from Earth, so using TT is justified:
  //
  class FlightTime
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    Time m_t;        // Time since Launch
    TT   m_launchTT; // Launch instant ("UnDef" if unknown or irrelevant)

  public:
    //-----------------------------------------------------------------------//
    // Ctors:                                                                //
    //-----------------------------------------------------------------------//
    // Only the Time Since Launch is given (if at all):
    //
    constexpr explicit FlightTime(Time a_t = Time())
    : m_t       (a_t),
      m_launchTT(TT::UnDef())
    {}

    // Both the Time Since Launch and the Launch Instant are given:
    //
    constexpr FlightTime(Time a_t, TT a_launch)
    : m_t       (a_t),
      m_launchTT(TT(a_launch))
    {}

    // Copy Ctor is auto-generated:
    constexpr FlightTime(FlightTime const&) = default;

    //-----------------------------------------------------------------------//
    // Assignment and Comparison:                                            //
    //-----------------------------------------------------------------------//
    // These ops are well-defined if both operands have the same  finite Launch
    // Instant, or one of the Launch Instants is UnDefined     (then the Launch
    // Instant of the result is UnDefined as well):
    //
    constexpr FlightTime& operator=  (FlightTime a_right)
    {
      bool   hasUnDef =  m_launchTT.IsUnDef() || a_right.m_launchTT.IsUnDef();
      assert(hasUnDef || m_launchTT == a_right.m_launchTT);
      m_t           = a_right.m_t;
      if (hasUnDef)
        m_launchTT  = TT::UnDef();
      return *this;
    }

    constexpr bool operator== (FlightTime a_right) const
    {
      // Check whether the args are comparable at all:
      CheckLaunchTTs(a_right);
      return m_t ==  a_right.m_t;
    }

    constexpr bool operator!= (FlightTime a_right) const
      { return !(*this == a_right); }

    constexpr bool operator>  (FlightTime a_right) const
    {
      // Check whether the args are comparable at all:
      CheckLaunchTTs(a_right);
      return m_t >   a_right.m_t;
    }

    constexpr bool operator>= (FlightTime a_right) const
    {
      // Check whether the args are comparable at all:
      CheckLaunchTTs(a_right);
      return m_t >=  a_right.m_t;
    }

    constexpr bool operator<  (FlightTime a_right) const
    {
      // Check whether the args are comparable at all:
      CheckLaunchTTs(a_right);
      return m_t <   a_right.m_t;
    }

    constexpr bool operator<= (FlightTime a_right) const
    {
      // Check whether the args are comparable at all:
      CheckLaunchTTs(a_right);
      return m_t <=  a_right.m_t;
    }

    //-----------------------------------------------------------------------//
    // Arithmetic:                                                           //
    //-----------------------------------------------------------------------//
    constexpr FlightTime& operator+= (Time a_dt)
    {
      m_t += a_dt;
      return *this;
    }

    constexpr FlightTime& operator-= (Time a_dt)
    {
      m_t -= a_dt;
      return *this;
    }

    constexpr  FlightTime  operator+ (Time a_dt) const
      { return FlightTime(m_t + a_dt, m_launchTT); }

    constexpr  FlightTime  operator- (Time a_dt) const
      { return FlightTime(m_t - a_dt, m_launchTT); }

    // Other way round:
    constexpr friend FlightTime operator+ (Time a_dt, FlightTime a_ft)
      { return FlightTime(a_ft.m_t + a_dt, a_ft.m_launchTT); }

    // Again, the difference of "FlightTime"s is only defined if both operands
    // have the same Launch Instant,  or if at least one of them has UnDefined
    // Launch Instant:
    constexpr Time operator- (FlightTime a_right) const
    {
      CheckLaunchTTs(a_right);
      return m_t - a_right.m_t;
    }

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    // Convertion to  "Time" yields Time Since Launch:
    constexpr operator Time  () const { return m_t;        }

    // Conversion to "TT" yields the absolute time; it requires LaunchTT to be
    // valid:
    constexpr operator TT    () const
    {
      assert(IsFinite(m_launchTT.GetTimeSinceEpoch()));
      return m_launchTT + m_t;
    }

    constexpr TT GetLaunchTT () const { return m_launchTT; }

  private:
    //-----------------------------------------------------------------------//
    // Utils:                                                                //
    //-----------------------------------------------------------------------//
    constexpr void CheckLaunchTTs([[maybe_unused]] FlightTime a_right) const
    {
      assert(m_launchTT.IsUnDef() || a_right.m_launchTT.IsUnDef() ||
             m_launchTT == a_right.m_launchTT);
    }
  };
};
// End namespace SpaceBallistics
