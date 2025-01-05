// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/CoOrds/TimeScales.hpp":                 //
//===========================================================================//
#pragma   once
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <utility>
#include <algorithm>

namespace SpaceBallistics
{
  namespace TBits
  {
    //-----------------------------------------------------------------------//
    // Table of Leap Seconds: [(Year, Month, Day)]:                          //
    //-----------------------------------------------------------------------//
    // NB: No more LeapSeconds are expected in the foreseeable future:
    //
    constexpr inline int  NLS = 27;
    constexpr inline std::tuple<int,int,int> LeapSecondDates[NLS]
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
  }

  //=========================================================================//
  // "UTC" Methods:                                                          //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "UTC" Non-Default Ctor (0):                                             //
  //-------------------------------------------------------------------------//
  constexpr UTC::UTC
  (
    std::tuple<int, int, int>    const& a_date,
    std::tuple<int, int, double> const& a_time
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

  //-------------------------------------------------------------------------//
  // "UTC" Non-Default Ctor (1):                                             //
  //-------------------------------------------------------------------------//
  constexpr UTC::UTC
  (
    int a_year, int a_month, int    a_day,
    int a_hour, int a_min,   double a_sec
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

  //-------------------------------------------------------------------------//
  // "UTC" Comparison:                                                       //
  //-------------------------------------------------------------------------//
  constexpr bool UTC::operator<  (UTC const& a_r) const
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

  //-------------------------------------------------------------------------//
  // "UTC::IsLeapYear":                                                      //
  //-------------------------------------------------------------------------//
  // Is this year a Gregorian Leap Year?
  //
  constexpr bool UTC::IsLeapYear(int a_year)
  {
    // It is assumed that the Gregorian Calendar is extended indefinitely into
    // the past (proleptic) and into the future across all jurisdictions. Yet
    // for technical reasons we need the Year to be positive:
    assert(a_year > 0);
    return a_year % 4 == 0 && (a_year % 100 != 0 || a_year % 400 == 0);
  }

  //-------------------------------------------------------------------------//
  // "UTC::DaysInMonth":                                                     //
  //-------------------------------------------------------------------------//
  constexpr int UTC::DaysInMonth(int a_year, int a_month)
  {
    assert(1 <= a_month && a_month <= 12);
    return
      (a_month == 4 || a_month == 6 || a_month == 9 || a_month == 11)
      ? 30 :
      (a_month == 2)
      ? (IsLeapYear(a_year) ? 29 : 28)
      : 31;
  }

  //-------------------------------------------------------------------------//
  // "UTC::HasLeapSecond":                                                   //
  //-------------------------------------------------------------------------//
  // Whether this "UTC" instant (considred up to a Minute) contains a Leap
  // Second:
  //
  constexpr bool UTC::HasLeapSecond() const
  {
    // A Leap Second can only occur at Hour=23, Min=59:
    std::tuple<int,int,int> date    { m_year, m_month, m_day };
    return m_hour == 23 && m_min == 59 &&
           (std::find
            (TBits::LeapSecondDates,
             TBits::LeapSecondDates + TBits::NLS, date) !=
             TBits::LeapSecondDates + TBits::NLS);
  }

  //-------------------------------------------------------------------------//
  // "UTC::IsValid": Over-All UTC Validation:                                //
  //-------------------------------------------------------------------------//
  constexpr bool UTC::IsValid() const
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

  //-------------------------------------------------------------------------//
  // "UTC::GetCumPrevLeapSecs":                                              //
  //-------------------------------------------------------------------------//
  // Cumulative number of Leap Seconds BEFORE this instant. If it corresponds
  // to a Leap Second itself, the current Leap Second is NOT counted yet:
  //
  constexpr Time UTC::GetCumPrevLeapSecs() const
  {
    int  res   = TBits::NLS;
    for (int i = TBits::NLS-1; i >= 0; --i, --res)
    {
      // Hear we use the fact that Leap Seconds always occur at the end of a
      // month:
      auto [y, m, _] =   TBits::LeapSecondDates[i];
      if (m_year > y || (m_year == y && m_month > m))
        break;
    }
    assert(res >= 0);
    return Time(double(res));
  }

  //=========================================================================//
  // Utils for UTC <-> TT Conversion:                                        //
  //=========================================================================//
  namespace TBits
  {
    // Fixed TT-TAI  diff:
    constexpr inline Time TT_TAI = 32.184_sec;

    //-----------------------------------------------------------------------//
    // Old-Style (1961--1972) "DeltaAT":                                     //
    //-----------------------------------------------------------------------//
    // DeltaAT = (TAI-UTC) difference, in sec.
    // The data are from the Explanatory Supplement to the Astronomical
    // Almanach, 3rd ed., 2013:
    //
    // "Old DataAT Node":
    struct ODATNode
    {
      UTC       m_until;
      Time_sec  m_a;
      Time_day  m_jd0; // JD_UTC
      Time_sec  m_b;
    };

    constexpr int      NODATN = 14;
    constexpr ODATNode ODATNodes[NODATN]
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
      Time_day  jd_utc(double(jd12h) - 0.5);

      // Apply the intra-day time:
      jd_utc += Time_day(double(a_utc.m_hour) /    24.0 +
                         double(a_utc.m_min)  /  1440.0 +
                         a_utc.m_sec          / 86400.0);

      // Finally, DeltaAT = TAI-UTC (initially 0):
      Time DeltaAT;

      // Try "Old-Style" nodes first:
      bool found = false;
      for (int i = 0; i < NODATN; ++i)
      {
        ODATNode const& node =  ODATNodes[i];
        // Check the monotonicity:
        assert(i == NODATN-1 || node.m_until < ODATNodes[i+1].m_until);

        if (a_utc < node.m_until)
        {
          DeltaAT =
            node.m_a + double((jd_utc - node.m_jd0) / 1.0_day) * node.m_b;
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
      return std::make_pair(jd_utc, DeltaAT);
    }

    //=======================================================================//
    // "MkMJSTT":                                                            //
    //=======================================================================//
    constexpr Time MkMJSTT(UTC const& a_utc)
    {
      // Compute the TAI components:
      auto [jd_utc, DeltaAT]  = MkTAI(a_utc);

      // Don't forget to include the fixed TT-TAI offset:
      return To_Time(jd_utc - Epoch_J2000) + DeltaAT + TT_TAI;
    }

    //=======================================================================//
    // Consts for TT -> UTC Conversion:                                      //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Old-Style (1961--1972) "DeltaAT" Nodes as MJS_TT:                     //
    //-----------------------------------------------------------------------//
    // XXX: a Boliler-Plate:
    //
    constexpr inline Time ODATNodesMJS[TBits::NODATN]
    {
      MkMJSTT(ODATNodes[ 0].m_until),
      MkMJSTT(ODATNodes[ 1].m_until),
      MkMJSTT(ODATNodes[ 2].m_until),
      MkMJSTT(ODATNodes[ 3].m_until),
      MkMJSTT(ODATNodes[ 4].m_until),
      MkMJSTT(ODATNodes[ 5].m_until),
      MkMJSTT(ODATNodes[ 6].m_until),
      MkMJSTT(ODATNodes[ 7].m_until),
      MkMJSTT(ODATNodes[ 8].m_until),
      MkMJSTT(ODATNodes[ 9].m_until),
      MkMJSTT(ODATNodes[10].m_until),
      MkMJSTT(ODATNodes[11].m_until),
      MkMJSTT(ODATNodes[12].m_until),
      MkMJSTT(ODATNodes[13].m_until)
    };

    //-----------------------------------------------------------------------//
    // Starts of all Leap Seconds as MJS_TT:                                 //
    //-----------------------------------------------------------------------//
    // XXX: a Boiler-Plate again:
    //
    constexpr inline std::tuple<int, int, double> LeapSecTime =
      std::make_tuple(23, 59, 60.0);

    constexpr inline Time LeapSecondsStartMJS[NLS]
    {
      MkMJSTT(UTC(LeapSecondDates[ 0], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 1], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 2], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 3], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 4], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 5], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 6], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 7], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 8], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[ 9], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[10], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[11], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[12], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[13], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[14], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[15], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[16], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[17], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[18], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[19], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[20], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[21], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[22], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[23], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[24], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[25], LeapSecTime)),
      MkMJSTT(UTC(LeapSecondDates[26], LeapSecTime))
    };
  }

  //=========================================================================//
  // Ctor: TT from UTC:                                                      //
  //=========================================================================//
  constexpr TT::TT(UTC const& a_utc)
  : m_MJS(TBits::MkMJSTT(a_utc))
  {}

  //=========================================================================//
  // TT -> UTC Conversion:                                                   //
  //=========================================================================//
  constexpr TT::operator UTC() const
  {
    //-----------------------------------------------------------------------//
    // Get the JD_UTC:                                                       //
    //-----------------------------------------------------------------------//
    Time_day jd_utc;  // Initially 0, which is never valid

    // Very old instants (BEFORE 1961): DeltaAT=0, so just subtract the fixed
    // TT_TAI offset:
    if (m_MJS < TBits::ODATNodesMJS[0])
      jd_utc = Epoch_J2000 + To_Time_day(m_MJS - TBits::TT_TAI);
    else
    {
      // Search "Old-Style" (1961--1972) Nodes:
      for (int i = 1; i < TBits::NODATN; ++i)   // NB: i==0 done above
      {
        // Check the monotonicity:
        assert(i == TBits::NODATN-1 ||
               TBits::ODATNodesMJS[i] < TBits::ODATNodesMJS[i+1]);

        // XXX: This "<" is susceptable to rounding errors!
        if (m_MJS < TBits::ODATNodesMJS[i])
        {
          TBits::ODATNode const& node = TBits::ODATNodes[i];
          // Again, check the monotonicity:
          assert(i == TBits::NODATN-1 ||
                 node.m_until < TBits::ODATNodes[i+1].m_until);

          // Solve the following equation wrt "jd_utc":
          // m_MJS   = To_Time(jd_utc - Epoch_J2000) + DeltaAT + TT_TAI
          // where
          // DeltaAT = node.m_a + node.m_b * (jd_utc - mode.m_jd0) / 1d :
          //
          Time denom = 86400.0_sec + node.m_b;
          double  c0 = double(86400.0_sec       / denom);
          double  c1 = double(node.m_b          / denom);
          double  c2 = double((m_MJS - node.m_a - TBits::TT_TAI)  / denom);
          jd_utc     = c0 * Epoch_J2000 + c1 * node.m_jd0 + c2 * 1.0_day;
          break;
        }
      }
    }
    // "jd_utc" has not been constructed? Then search the Leap Seconds:
    //
    bool isLS  = false; // m_MJS is in [LeapSecBegin .. LeapSecEnd) ?
    if (IsZero(jd_utc))
    {
      // Then we are from 1972 on, so Leap Seconds are to be used:
      assert(m_MJS >= TBits::ODATNodesMJS[TBits::NODATN-1]);

      // Search the "TBits::LeapSecondsStartMJS" table:
      for (int i = 0; i < TBits::NLS; ++i)
      {
        // End of the "i"th Leap Second:
        Time startI = TBits::LeapSecondsStartMJS[i];
        Time endI   = startI + 1.0_sec;

        // Check the monotonicity:
        assert(i == TBits::NLS-1 || startI < TBits::LeapSecondsStartMJS[i+1]);

        // XXX: This "<" is susceptable to rounding errors!
        if (m_MJS < endI)
        {
          // Previous COMPLETE Leap Seconds (i) and the fixed offset:
          Time off = Time(double(i)) + 10.0_sec + TBits::TT_TAI;

          // Then "jd_utc" is:
          jd_utc = Epoch_J2000 + To_Time_day(m_MJS - off);

          // Are we are actually within the "i"th Leap Second (incomplete
          // yet)?
          isLS = (startI <= m_MJS);
          break;
        }
      }
      // "jd_utc" not constructed yet?
      if (IsZero(jd_utc))
      {
        // Then all currently-known Leap Seconds have been traversed, and are
        // in the past. Then the last available Number of Leap Seconds is ext-
        // rapolated:
        assert(m_MJS >= TBits::LeapSecondsStartMJS[TBits::NLS-1] + 1.0_sec);
        Time   off    = Time(double(TBits::NLS)) + 10.0_sec + TBits::TT_TAI;
        jd_utc        = Epoch_J2000 + To_Time_day(m_MJS - off);
      }
    }
    // We must have obtained valid "jd_utc" and possibly the curr Leap Second
    // info. Ante-Deluvian dates are not allowed:
    assert(IsPos(jd_utc));

    // Extract the Gregorian Calendar Date (possibly proleptic) from "jd_utc":
    // From:
    // Fliegel, H. & Van Flandern, T.  Comm. of the ACM, Vol. 11, No. 10,
    // October 1968, p. 657:
    //
    jd_utc         += 0.5_day;
    // XXX: "Floor" is susceptable to rounding errors!
    double   jdW    = Floor(double(jd_utc / 1.0_day));
    long     jd     = long (jdW);
    Time     daySec = To_Time(jd_utc - Time_day(jdW));
    assert(!IsNeg(daySec));

    long k     = jd + 68569;
    long n     = (4 * k) / 146097;
    k         -= (146097 *  n + 3)  / 4;
    long m     = (  4000 * (k + 1)) / 1461001;
    k         -= (1461 * m) / 4 - 31;

    // Year, Month, Day:
    long lm    = (80 * k) / 2447;
    int day    = int(k - (2447 * lm) / 80);
    k          = lm / 11;
    int month  = int(lm  + 2  - 12  * k);
    int year   = int(100 * (n - 49) + m + k);
    assert(year > 0 && 1 <= month && month <= 12 && 1 <= day &&
           day <= UTC::DaysInMonth(year, month));

    // Hour, Min, Sec:
    // XXX: Again,  "Floor" is susceptable to rounding errors!
    int   iSec = int(Floor(double(daySec / 1.0_sec)));
    int   hour = iSec / 3600;
    int    min = iSec / 60 - hour * 60;
    double sec = double
                 ((daySec - Time(double(3600 * hour + 60 * min))) / 1.0_sec);
    assert(0 <= hour  && hour <= 23 && 0 <= min && min <= 59 &&
           0.0 <= sec && sec  <  60.0);

    // If it is a LeapSecond, it would as yet be 0th of the following minute,
    // but we need to mark it as the 60th, and adjust the minute, hour and
    // day:
    if (isLS)
    {
      assert(Floor(sec) == 0.0 && min == 0 && hour == 0 && day == 1 &&
            (month == 1 || month == 7));
      sec   = 60.0 + (sec - Floor(sec));
      min   = 59;
      hour  = 23;
      if (month == 7)
      {
        day   = 30;
        month = 6;
        // year unchanged
      }
      else
      {
        day   = 31;
        month = 12;
        --year;
      }
    }
    // We can now compose the UTC obj:
    return UTC(year, month, day, hour, min, sec);
  }
};
// End namespace SpaceBallistics
