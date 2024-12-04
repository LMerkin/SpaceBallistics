// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/TimeScales.h":                  //
//===========================================================================//
#pragma   once
#include "SpaceBallistics/Types.hpp"
#include <algorithm>
#include <tuple>
#include <utility>

namespace SpaceBallistics
{
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

    //-----------------------------------------------------------------------//
    // Consts:                                                               //
    //-----------------------------------------------------------------------//
    // Table of Leap Seconds: [(Year, Month, Day)]. NB: No more LeapSeconds are
    // expected in the foreseeable future:
    constexpr static int  NLS = 27;
    constexpr static std::tuple<int,int,int> LeapSecondDates[NLS]
    {
      { 1972,  6, 30 }, { 1972, 12, 31 }, { 1973, 12, 31 },
      { 1974, 12, 31 }, { 1975, 12, 31 }, { 1976, 12, 31 },
      { 1977, 12, 31 }, { 1978, 12, 31 }, { 1979, 12, 31 },
      { 1981,  6, 30 }, { 1982,  6, 30 }, { 1983,  6, 30 },
      { 1985,  6, 30 }, { 1987, 12, 31 }, { 1989, 12, 31 },
      { 1990, 12, 31 }, { 1992,  6, 30 }, { 1993,  6, 30 },
      { 1994,  6, 30 }, { 1995, 12, 31 }, { 1997,  6, 30 },
      { 1998, 12, 31 }, { 2005, 12, 31 }, { 2008, 12, 31 },
      { 2012,  6, 30 }, { 2015,  6, 30 }, { 2016, 12, 31 }
    };

    //-----------------------------------------------------------------------//
    // "IsLeapYear":                                                         //
    //-----------------------------------------------------------------------//
    // Is this year a Gregorian Leap Year?
    //
    constexpr bool IsLeapYear() const
      { return  m_year % 4 == 0 && (m_year % 100 != 0 || m_year % 400 == 0); }

    //-----------------------------------------------------------------------//
    // "HasLeapSecond":                                                      //
    //-----------------------------------------------------------------------//
    // Whether this (Year, Month, Day, Hour, Min) contains a Leap Second:
    //
    constexpr bool HasLeapSecond() const
    {
      // A Leap Second can only occur at Hour=23, Min=59:
      std::tuple<int,int,int> date    { m_year, m_month, m_day };
      return m_hour == 23 && m_min == 59 &&
             std::find(LeapSecondDates, LeapSecondDates + NLS, date) !=
             LeapSecondDates + NLS;
    }

    //-----------------------------------------------------------------------//
    // "IsValid": Over-All UTC Validation:                                   //
    //-----------------------------------------------------------------------//
    constexpr bool IsValid() const
    {
      int daysInMonth =
          (m_month == 4 || m_month == 6 || m_month == 9 || m_month == 11)
          ? 30 :
          (m_month == 2)
          ? (IsLeapYear() ? 29 : 28)
          : 31;

      double secsInMin = HasLeapSecond() ? 61.0 : 60.0;

      // NB: We currently require m_year > 0:
      return
        0   <  m_year                            &&
        1   <= m_month && m_month <= 12          &&
        1   <= m_day   && m_day   <= daysInMonth &&
        0   <= m_hour  && m_hour  <= 23          &&
        0   <= m_min   && m_min   <= 59          &&
        0.0 <= m_sec   && m_sec   <  secsInMin;
    }

    //-----------------------------------------------------------------------//
    // "CumPrevLeapSecs":                                                    //
    //-----------------------------------------------------------------------//
    // Cumulative number of Leap Seconds BEFORE this instant. If it corresponds
    // to a Leap Second itself, the current Leap Second does NOT count yet:
    //
    constexpr Time CumPrevLeapSecs() const
    {
      int res = NLS;
      for (int i = NLS-1; i >= 0; --i, --res)
      {
        // Hear we use the fact that Leap Seconds always occur at the end of a
        // month:
        auto [y, m, _] =  LeapSecondDates[i];
        if (m_year > y || (m_year == y && m_month > m))
          break;
      }
      assert(res >= 0);
      return Time(double(res));
    }
  };

  //=========================================================================//
  // "TT" Class:                                                             //
  //=========================================================================//
  // Terrestrial Time (Uniform, Relativistic), implemented via TAI (with an
  // offset):
  class TDB;

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
    // NB: Default Ctor returns the Epoch (J2000.0):
    constexpr TT              ()                = default;
    constexpr TT              (TT const&)       = default;
    constexpr TT&  operator=  (TT const&)       = default;
    constexpr bool operator== (TT const&) const = default;

    // The J2000.0 Epoch = 2000 Jan 1.5, in TT:
    constexpr static Time_day Epoch    = 2451545.0_day;
    constexpr static Time_jyr EpochJYr = Time_jyr(2000.0 + 0.5 / 365.25);

    // Fixed TT-TAI offset:
    constexpr static Time     TT_TAI = 32.184_sec;

    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // From UTC:
    constexpr explicit TT(UTC const& a_utc)
    {
      // Compute the TAI components:
      auto [JD_UTC, DeltaAT]  = MkTAI(a_utc);

      // Don't forget to include the fixed TT-TAI offset:
      m_MJS  = To_Time(JD_UTC - Epoch) + DeltaAT + TT_TAI;
    }

    // Directly constructing TT from JD_TT:
    constexpr explicit TT(Time_day a_jd_tt)
    : m_MJS    (To_Time_sec(a_jd_tt - Epoch))
    {}

    // Directly constructing from JYear_TT. XXX: In this case we assume that
    // the Epoch is 2000.0 exactly. XXX: This is for "low-precision" computa-
    // tions such as Precession Matrices:
    //
    constexpr explicit TT(Time_jyr a_jyr_tt)
    : m_MJS    (To_Time_sec( a_jyr_tt - EpochJYr))
    {}

    // From TDB (not "constexpr", as DE440T is required):
    explicit TT(TDB const& a_tdb);

    //-----------------------------------------------------------------------//
    // Extracting the Time value (since the Epoch) or JD_TT val:             //
    //-----------------------------------------------------------------------//
    // XXX: USE WITH CARE, as TimeScale info is then lost:
    //
    constexpr Time     GetTimeSinceEpoch() const
      { return m_MJS; }

    constexpr Time_day GetJDsSinceEpoch () const
      { return To_Time_day(m_MJS); }

    // Extracting the JD_TDB:
    constexpr Time_day GetJD            () const
      { return Epoch + To_Time_day(m_MJS); }

    //-----------------------------------------------------------------------//
    // "MkTAI":                                                              //
    //-----------------------------------------------------------------------//
    // Util used as part of UTC->TT conversion. For a given UTC, returns the
    // pair (JD_UTC, DeltaAT), where DeltaAT = TAI - UTC. Can be used in the
    // stand-alone mode as well:
    //
    constexpr static std::pair<Time_day, Time> MkTAI(UTC const& a_utc)
    {
      // Verify the arg:
      assert(a_utc.IsValid());

      // First, compute JD_UTC for the Date part only, including Gregorian Leap
      // Leap Years (retroactively applied if necessary).
      // The algorithm is from NOVAS-C Ver 3.1, originally by:
      // Fliegel, H. & Van Flandern, T., Comm. of the ACM, Vol.11, No.10,
      // October 1968, p.657:
      //
      int mm    = (a_utc.m_month - 14) / 12;

      // This assumes that JD starts @ 12 noon:
      int jd12h =
        a_utc.m_day
        - 32075
        + (1461 * ( a_utc.m_year  + 4800 + mm))        / 4
        + ( 367 * ( a_utc.m_month -    2 - mm  *  12)) / 12
        - (   3 * ((a_utc.m_year  + 4900 + mm) / 100)) / 4;

      // Convert the JD to mid-night start:
      Time_day  JD_UTC(double(jd12h) - 0.5);

      // Apply the intra-day time:
      JD_UTC += Time_day(double(a_utc.m_hour) /    24.0 +
                         double(a_utc.m_min)  /  1440.0 +
                         a_utc.m_sec          / 86400.0);

      // MJD (Modified Julian Days) is used in the DeltaAT computation below:
      Time_day mjd = JD_UTC - 2'400'000.5_day;

      // DeltaAT = (TAI-UTC) difference, in sec.
      // The formulas are from the Explanatory Supplement to the Astronomical
      // Almanach, 3rd ed., 2013:
      // XXX: We have to use "Magnitude" here:
      //
      Time DeltaAT =
        //
        // Prior to 1961: The diff is not known, XXX: assume it to be 0 :
        //
        (a_utc.m_year <  1961)
        ? 0.0_sec :
        //
        // Prior to introduction of Leap Seconds in 1972:
        //
        (a_utc.m_year == 1961 && a_utc.m_month < 8)
        ? 1.422818_sec + (mjd - 37300.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year <  1962)
        ? 1.372818_sec + (mjd - 37300.0_day).Magnitude() * 0.001296_sec 
        :
        (a_utc.m_year == 1962 || (a_utc.m_year == 1963 && a_utc.m_month < 11))
        ? 1.845858_sec + (mjd - 37665.0_day).Magnitude() * 0.0011232_sec
        :
        (a_utc.m_year == 1963)
        ? 1.945858_sec + (mjd - 37665.0_day).Magnitude() * 0.0011232_sec
        :
        (a_utc.m_year == 1964 && a_utc.m_month < 4)
        ? 3.240130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1964 && a_utc.m_month < 9)
        ? 3.340130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1964)
        ? 3.440130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1965 && a_utc.m_month < 3)
        ? 3.540130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1965 && a_utc.m_month < 7)
        ? 3.640130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1965 && a_utc.m_month < 9)
        ? 3.740130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year == 1965)
        ? 3.840130_sec + (mjd - 38761.0_day).Magnitude() * 0.001296_sec
        :
        (a_utc.m_year <  1968 || (a_utc.m_year == 1968 && a_utc.m_month < 2))
        ? 4.313170_sec + (mjd - 39126.0_day).Magnitude() * 0.002592_sec
        :
        (a_utc.m_year <  1972)
        ? 4.213170_sec + (mjd - 39126.0_day).Magnitude() * 0.002592_sec
        :
        // Since 1972, the Leap Seconds are used (don't forget the initial
        // offset of 10 sec!):
        10.0_sec + a_utc.CumPrevLeapSecs();

        // All Done:
        return std::make_pair(JD_UTC, DeltaAT);
    }

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
    // NB: Default Ctor returns the Epoch (J2000.0):
    constexpr TDB              ()                 = default;
    constexpr TDB              (TDB const&)       = default;
    constexpr TDB&  operator=  (TDB const&)       = default;
    constexpr bool  operator== (TDB const&) const = default;

    // The J2000.0 Epoch = 2000 Jan 1.5, in TDB:
    constexpr static Time_day Epoch = 2451545.0_day;
    
    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // Constructing TDB from TT. The implementation is non-trivial, requires
    // the JPL Ephemerides (DE440T), so it is not a "constexpr":
    explicit TDB(TT const& a_tt);

    // Directly constructing TDB from JD_TDB:
    constexpr explicit TDB (Time_day a_jd_tdb)
    : m_MJS    (To_Time_sec(a_jd_tdb - Epoch))
    {}

    //-----------------------------------------------------------------------//
    // Extracting the Time value (since the Epoch) or JD_TDB val:            //
    //-----------------------------------------------------------------------//
    // XXX: USE WITH CARE, as TimeScale info is then lost:
    //
    constexpr Time     GetTimeSinceEpoch() const
      { return m_MJS; }

    constexpr Time_day GetJDsSinceEpoch () const
      { return To_Time_day(m_MJS); }

    // Extracting the JD_TDB:
    constexpr Time_day GetJD            () const
      { return Epoch + To_Time_day(m_MJS); }

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
  };
};
// End namespace SpaceBallistics
