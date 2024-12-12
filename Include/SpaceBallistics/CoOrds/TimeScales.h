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

    // NB: Default Ctor is implicitly deleted!

    //-----------------------------------------------------------------------//
    // Non-Default Ctors: From Date and Time:                                //
    //-----------------------------------------------------------------------//
    constexpr explicit UTC
    (
      std::tuple<int, int, int>    const& a_date,
      std::tuple<int, int, double> const& a_time = std::make_tuple(0, 0, 0.0)
    )
    : m_year (std::get<0>(a_date)),
      m_month(std::get<1>(a_date)),
      m_day  (std::get<2>(a_date)),
      m_hour (std::get<0>(a_time)),
      m_min  (std::get<1>(a_time)),
      m_sec  (std::get<2>(a_time))
    {
      // Verification:
      assert(IsValid());
    }

    constexpr explicit UTC
    (
      int a_year,     int a_month = 1, int    a_day = 1,
      int a_hour = 0, int a_min   = 0, double a_sec = 0.0
    )
    : m_year (a_year),
      m_month(a_month),
      m_day  (a_day),
      m_hour (a_hour),
      m_min  (a_min),
      m_sec  (a_sec)
    {
      // Verification:
      assert(IsValid());
    }

    //-----------------------------------------------------------------------//
    // Comparison:                                                           //
    //-----------------------------------------------------------------------//
    constexpr bool operator== (UTC const& a_r) const = default;

    constexpr bool operator<  (UTC const& a_r) const
    {
      return
        (m_year <  a_r.m_year)
        ||
        (m_year == a_r.m_year && m_month <  a_r.m_month)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  <  a_r.m_day)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  <  a_r.m_day)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  == a_r.m_day  && m_hour  <  a_r.m_hour)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  == a_r.m_day  && m_hour  <  a_r.m_hour)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  == a_r.m_day  && m_hour  == a_r.m_hour   &&
         m_min  <  a_r.m_min)
        ||
        (m_year == a_r.m_year && m_month == a_r.m_month  &&
         m_day  == a_r.m_day  && m_hour  == a_r.m_hour   &&
         m_min  == a_r.m_min  && m_sec   <  a_r.m_sec);
    }

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
    constexpr static bool IsLeapYear(int a_year)
    {
      // It is assumed that the Gregorian Calendar is extended indefinitely into
      // the past (and into the future) across all jurisdictions. Yet for techn-
      // ical reasons we need the Year to be positive:
      assert(a_year > 0);
      return a_year % 4 == 0 && (a_year % 100 != 0 || a_year % 400 == 0);
    }

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
    // "DaysInMonth":                                                        //
    //-----------------------------------------------------------------------//
    constexpr static int DaysInMonth(int a_year, int a_month)
    {
      assert(1 <= a_month && a_month <= 12);
      return
        (a_month == 4 || a_month == 6 || a_month == 9 || a_month == 11)
        ? 30 :
        (a_month == 2)
        ? (IsLeapYear(a_year) ? 29 : 28)
        : 31;
    }

    //-----------------------------------------------------------------------//
    // "IsValid": Over-All UTC Validation:                                   //
    //-----------------------------------------------------------------------//
    constexpr bool IsValid() const
    {
      double secsInMin = HasLeapSecond() ? 61.0 : 60.0;

      // NB: We currently require m_year > 0:
      return
        0   <  m_year                                             &&
        1   <= m_month && m_month <= 12                           &&
        1   <= m_day   && m_day   <= DaysInMonth(m_year, m_month) &&
        0   <= m_hour  && m_hour  <= 23                           &&
        0   <= m_min   && m_min   <= 59                           &&
        0.0 <= m_sec   && m_sec   <  secsInMin;
    }

    //-----------------------------------------------------------------------//
    // "GetCumPrevLeapSecs":                                                 //
    //-----------------------------------------------------------------------//
    // Cumulative number of Leap Seconds BEFORE this instant. If it corresponds
    // to a Leap Second itself, the current Leap Second is NOT counted yet:
    //
    constexpr Time GetCumPrevLeapSecs() const
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
  // Utils for UTC <-> TT Conversion:                                        //
  //=========================================================================//
  namespace TBits
  {
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // The J2000.0 Epoch = 2000 Jan 1.5 (will be used in TT and TDB):
    //
    constexpr inline Time_day Epoch_J2000    = 2451545.0_day;
    constexpr inline Time_jyr Epoch_J2000_Yr = Time_jyr(2000.0 + 0.5 / 365.25);

    // Fixed TT-TAI  diff:
    constexpr inline Time     TT_TAI         = 32.184_sec;

    //=======================================================================//
    // Old-Style (1961--1972) "DeltaAT" Nodes:                               //
    //=======================================================================//
    struct OldDeltaATNode
    {
      UTC       m_utc;
      Time_sec  m_a;
      Time_day  m_base;
      Time_sec  m_b;
    };

    // DeltaAT = (TAI-UTC) difference, in sec.
    // The data are from the Explanatory Supplement to the Astronomical
    // Almanach, 3rd ed., 2013:
    constexpr int            NODAT   = 14;
    constexpr OldDeltaATNode ODATNodes[NODAT]
    {
      // Prior to 1961: The diff is not known, XXX: assume it to be 0:
      { UTC(1961,  1), 0.0_sec,      2437300.5_day, 0.0_sec       },
      { UTC(1961,  8), 1.422818_sec, 2437300.5_day, 0.001296_sec  },
      { UTC(1962,  1), 1.372818_sec, 2437300.5_day, 0.001296_sec  },
      { UTC(1963, 11), 1.845858_sec, 2437665.5_day, 0.0011232_sec },
      { UTC(1964,  1), 1.945858_sec, 2437665.5_day, 0.0011232_sec },
      { UTC(1964,  4), 3.240130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1964,  9), 3.340130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1965,  1), 3.440130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1965,  3), 3.540130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1965,  7), 3.640130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1965,  9), 3.740130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1966,  1), 3.840130_sec, 2438761.5_day, 0.001296_sec  },
      { UTC(1968,  2), 4.313170_sec, 2439126.5_day, 0.002592_sec  },
      { UTC(1972,  1), 4.213170_sec, 2439126.5_day, 0.002592_sec  }
      // Since 1972, Leap Seconds are used...
    };

    //=======================================================================//
    // "MkTAI":                                                              //
    //=======================================================================//
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

      // Finally, DeltaAT = TAI-UTC (initially 0):
      Time DeltaAT;

      // Try "Old-Style" nodes first:
      bool found = false;
      for (int i = 0; i < NODAT; ++i)
      {
        OldDeltaATNode const& node = ODATNodes[i];
        if (a_utc < node.m_utc)
        {
          DeltaAT =
            node.m_a + double((JD_UTC - node.m_base) / 1.0_day) * node.m_b;
          found   = true;
          break;
        }
      }
      // If we got here and "DeltaAT" has not neen found yet,
      // so UTC must be >= 1972: Use Leap Seconds:
      if (!found)
      {
        assert(a_utc.m_year >= 1972);
        // Don't forget the initial 10.0 Leap Seconds:
        DeltaAT = 10.0_sec + a_utc.GetCumPrevLeapSecs();
      }
      // All Done:
      return std::make_pair(JD_UTC, DeltaAT);
    }

    //=======================================================================//
    // "MkMJSTT":                                                            //
    //=======================================================================//
    constexpr Time MkMJSTT(UTC const& a_utc)
    {
      // Compute the TAI components:
      auto [JD_UTC, DeltaAT]  = MkTAI(a_utc);

      // Don't forget to include the fixed TT-TAI offset:
      return To_Time(JD_UTC - Epoch_J2000) + DeltaAT + TT_TAI;
    }
  }

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
    // NB: Default Ctor actually returns the Epoch_J2000 in TT:
    constexpr TT              ()                = default;
    constexpr TT              (TT const&)       = default;
    constexpr TT&  operator=  (TT const&)       = default;
    constexpr bool operator== (TT const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // From UTC:
    constexpr explicit  TT(UTC const& a_utc)
    : m_MJS(TBits::MkMJSTT(a_utc))
    {}

    // Directly constructing TT from JD_TT:
    constexpr explicit TT(Time_day a_jd_tt)
    : m_MJS    (To_Time_sec(a_jd_tt - TBits::Epoch_J2000))
    {}

    // Directly constructing from JYear_TT:
    constexpr explicit TT(Time_jyr a_jyr_tt)
    : m_MJS    (To_Time_sec( a_jyr_tt - TBits::Epoch_J2000_Yr))
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
      { return TBits::Epoch_J2000 + To_Time_day(m_MJS); }

    //-----------------------------------------------------------------------//
    // Consts for "GetUTC":                                                  //
    //-----------------------------------------------------------------------//
    // Old-Style (1961--1972) "DeltaAT" Nodes in MJS_TT:
    //
    constexpr static Time OldDeltaATNodesMJS[TBits::NODAT]
    {
      TBits::MkMJSTT(TBits::ODATNodes[ 0].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 1].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 2].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 3].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 4].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 5].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 6].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 7].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 8].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[ 9].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[10].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[11].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[12].m_utc),
      TBits::MkMJSTT(TBits::ODATNodes[13].m_utc)
    };

    // Beginnings of Leap Seconds in MJS_TT:
    //
    constexpr static std::tuple<int, int, double> LeapSecTime =
      std::make_tuple(23, 59, 60.0);

    constexpr static Time LeapSecondsStartMJS[UTC::NLS]
    {
      // XXX: A horrible boiler-plate, but there is no easy way to constract a
      // "constexpr" array via a mapping function. Fortunately, no new Leap Se-
      // conds are expected any time soon:
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 0], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 1], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 2], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 3], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 4], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 5], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 6], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 7], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 8], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[ 9], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[10], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[11], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[12], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[13], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[14], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[15], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[16], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[17], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[18], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[19], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[20], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[21], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[22], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[23], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[24], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[25], LeapSecTime)),
      TBits::MkMJSTT(UTC(UTC::LeapSecondDates[26], LeapSecTime))
    };

    //-----------------------------------------------------------------------//
    // "GetDeltaAT":                                                         //
    //-----------------------------------------------------------------------//
    // Unlike "TBits::GetDeltaAT", here the arg is TT:
    //
/*
    constexpr Time GetDeltaAT() const
    {
    }

    //-----------------------------------------------------------------------//
    // "GetUTC": TT -> UTC Conversion:                                       //
    //-----------------------------------------------------------------------//
    constexpr operator UTC() const
    {
      // First, try Old-Style Nodes:
    }
*/

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
    : m_MJS    (To_Time_sec(a_jd_tdb - TBits::Epoch_J2000))
    {}

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
      { return TBits::Epoch_J2000 + To_Time_day(m_MJS); }

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
