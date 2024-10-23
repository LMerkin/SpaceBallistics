// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/TimeScales.h":                  //
//===========================================================================//
#pragma   once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TT" Class:                                                             //
  //=========================================================================//
  // Terrestrial Time (Uniform, Relativistic), implemented via TAI:
  //
  class TT
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // The time is stored in "Julian Seconds", ie as a number of TAI seconds
    // since the Epoch:
    Time  m_JS;

  public:
    //-----------------------------------------------------------------------//
    // Default Ctor, Copy Ctor, Assignment, Equality:                        //
    //-----------------------------------------------------------------------//
    // Default Ctor is auto-generated though it makes little sense (corresponds
    // to JS=0):
    constexpr TT              ()                = default;
    constexpr TT              (TT const&)       = default;
    constexpr TT&  operator=  (TT const&)       = default;
    constexpr bool operator== (TT const&) const = default;

    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // Constructing TT from UTC:
    constexpr TT
    (
      int     a_year,     //  Currently 1650 .. 2150 (for DE440T compatibility)
      int     a_month,    //  1..12
      int     a_day,      //  1..DaysInMonth(Year, Month)
      int     a_hour,     //  0..23
      int     a_min,      //  0..59
      double  a_sec);     // [0..60) or [0..61) for Leap Seconds

    // Directly constructing TT from JD or Seconds since the Epoch:
    // XXX: USE WITH CARE!
    constexpr TT  (Time_day a_jd_tt)
    : m_JS     (To_Time_sec(a_jd_tt))
    {}

    constexpr TT  (Time a_js_tt)
    : m_JS(a_js_tt)
    {}

    // Extracting the Time value (since the Epoch): XXX: USE WITH CARE!
    constexpr Time GetTime() const { return m_JS; }

    //-----------------------------------------------------------------------//
    // Time Intervals (Durations):                                           //
    //-----------------------------------------------------------------------//
    constexpr Time operator-  (TT   a_right) const
      { return m_JS - a_right.m_JS;    }

    constexpr TT&  operator+= (Time a_right)
    {
      m_JS += a_right;
      return  *this;
    }

    constexpr TT&  operator-= (Time a_right)
    {
      m_JS -= a_right;
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
      { return m_JS <  a_right.m_JS; }

    constexpr bool operator<= (TT a_right) const
      { return m_JS <= a_right.m_JS; }

    constexpr bool operator>  (TT a_right) const
      { return m_JS >  a_right.m_JS; }

    constexpr bool operator>= (TT a_right) const
      { return m_JS >= a_right.m_JS; }
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
    // since the Epoch:
    Time  m_JS;

  public:
    //-----------------------------------------------------------------------//
    // Default Ctor, Copy Ctor, Assignment, Equality:                        //
    //-----------------------------------------------------------------------//
    // Default Ctor is auto-generated though it makes little sense (corresponds
    // to JS=0):
    constexpr TDB              ()                 = default;
    constexpr TDB              (TDB const&)       = default;
    constexpr TDB&  operator=  (TDB const&)       = default;
    constexpr bool  operator== (TDB const&) const = default;
    
    //-----------------------------------------------------------------------//
    // Non-Default Ctors:                                                    //
    //-----------------------------------------------------------------------//
    // Constructing TDB from TT. The implementation is non-trivial, requires
    // JPL Ephemerides (DE440T):
    //
    constexpr TDB(TT a_tt);

    // Directly constructing TDB from JD or Seconds since the Epoch:
    // XXX: USE WITH CARE!
    constexpr TDB (Time_day a_jd_tdb)
    : m_JS     (To_Time_sec(a_jd_tdb))
    {}

    constexpr TDB (Time a_js_tdb)
    : m_JS(a_js_tdb)
    {}

    // Extracting the Time value (since the Epoch): XXX: USE WITH CARE!
    constexpr Time GetTime() const { return m_JS; }

    //-----------------------------------------------------------------------//
    // Time Intervals (Durations):                                           //
    //-----------------------------------------------------------------------//
    constexpr Time operator-  (TDB  a_right) const
      { return m_JS - a_right.m_JS;    }
    
    constexpr TDB& operator+= (Time a_right)
    {
      m_JS += a_right;
      return  *this;
    }

    constexpr TDB& operator-= (Time a_right)
    {
      m_JS -= a_right;
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
      { return m_JS <  a_right.m_JS; }

    constexpr bool operator<= (TDB a_right) const
      { return m_JS <= a_right.m_JS; }

    constexpr bool operator>  (TDB a_right) const
      { return m_JS >  a_right.m_JS; }

    constexpr bool operator>= (TDB a_right) const
      { return m_JS >= a_right.m_JS; }
  };
};
// End namespace SpaceBallistics
