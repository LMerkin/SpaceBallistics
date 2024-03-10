// vim:ts=2:et
//===========================================================================//
//                       "SpaceBallistics/Utils.hpp":                        //
//                                Misc Utils                                 //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include <tuple>

namespace SpaceBallistics
{
  //=========================================================================//
  // "ToRads":                                                               //
  //=========================================================================//
  // Conversion of an angle given as ('+'|'-', Degs, Mins, Secs) into Radians:
  //
  constexpr inline double ToRads(std::tuple<char,int,int,double> const& a_angle)
  {
    char   sign = std::get<0>(a_angle);
    assert(sign == '+' || sign == ' ' || sign == '-');

    double degs = double(std::get<1>(a_angle));   // Any value allowed
    double mins = double(std::get<2>(a_angle));   // Must be 0..59
    assert(0.0 <= mins && mins <= 59.0);
    double secs = std::get<3> (a_angle);          // Must be 0..60-
    assert(0.0 <= secs && secs <  60.0);

    double rads =
      double(To_Angle_rad(Angle_deg (double(degs))) +
             To_Angle_rad(Angle_amin(double(mins))) +
             To_Angle_rad(Angle_asec(secs)));

    return (sign == '-') ? (-rads) : rads;
  }
}
