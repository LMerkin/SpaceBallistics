// vim:ts=2:et
//===========================================================================//
//                           "Tests/TestUtils.hpp":                          //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/TimeScales.hpp"
#include <cassert>
#include <ctime>

namespace
{
	using namespace SpaceBallistics;

  //=========================================================================//
  // "TTofStr":                                                              //
  //=========================================================================//
  inline TT TTofStr(char const* a_str)
  {
    assert(a_str != nullptr && *a_str != '\0');

    // First of all, if "a_str" contains the "JD" prefix, we assume it is TT
    // in the JD format:
    if (a_str[0] == 'J' && a_str[1] == 'D')
      return TT(Time_day(atof(a_str + 2)));

    // Otherwise, it must be a UTC time-stamp in the format
    // "YYYY-MM-DD_hh:mm:ss":
    tm tmUTC;
    if (*strptime(a_str, "%Y-%m-%d_%H:%M:%S", &tmUTC) != '\0')
      return TT();

    UTC utc
    (
      tmUTC.tm_year + 1900,
      tmUTC.tm_mon  + 1,
      tmUTC.tm_mday,
      tmUTC.tm_hour,
      tmUTC.tm_min,
      double(tmUTC.tm_sec)
    );
    return TT(utc);
  }

  //=========================================================================//
  // "TimeOfStr":                                                            //
  //=========================================================================//
  inline Time TimeOfStr(char const* a_str)
  {
    assert (a_str != nullptr && *a_str != '\0');
    switch (a_str[0])
    {
      case 's': return Time            (atof(a_str + 1));
      case 'd': return To_Time(Time_day(atof(a_str + 1)));
      case 'y': return To_Time(Time_jyr(atof(a_str + 1)));
      default : return Time();
    }
  }
}
