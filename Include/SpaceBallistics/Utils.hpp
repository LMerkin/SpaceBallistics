// vim:ts=2:et
//===========================================================================//
//                       "SpaceBallistics/Utils.hpp":                        //
//                                Misc Utils                                 //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include <tuple>

namespace DimTypes
{
  //=========================================================================//
  // Conversion of Angles:                                                   //
  //=========================================================================//
  namespace SB = SpaceBallistics;

  // Allowing explicit conversions of "Angle[_rad]" into "double":
  template<>
  constexpr inline
  DimQ
  <
    SB::DimQ_Encs::DimExp(unsigned(SB::DimsE::Angle)),
    SB::DimQ_Encs::MkUnit(unsigned(SB::DimsE::Angle), 0), // 0=rad: fund unit
    SB::DimQ_RepT,
    SB::DimQ_MaxDims
  >
  ::operator double() const
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle>   &&
                  std::is_same_v<ThisDimQ, SB::Angle_rad>);
    return double(*this / SB::Angle(1.0));
  }

  // Similarly, converting "Angle_deg" into "double", via "Angle" (rad):
  template<>
  constexpr inline
  DimQ
  <
    SB::DimQ_Encs::DimExp(unsigned(SB::DimsE::Angle)),
    SB::DimQ_Encs::MkUnit(unsigned(SB::DimsE::Angle), 1), // 1=deg
    SB::DimQ_RepT,
    SB::DimQ_MaxDims
  >
  ::operator double() const
  {
    using ThisDimQ = std::remove_cv_t<std::remove_reference_t<decltype(*this)>>;
    static_assert(std::is_same_v<ThisDimQ, SB::Angle_deg>);
    return double(SB::To_Angle(*this));
  }
}

namespace SpaceBallistics
{
  //=========================================================================//
  // "ToRads":                                                               //
  //=========================================================================//
  // Conversion of an angle given as ('+'|'-', Degs, Mins, Secs) into Radians:
  //
  constexpr inline Angle ToRads(std::tuple<char,int,int,double> const& a_angle)
  {
    char   sign = std::get<0>(a_angle);
    assert(sign == '+' || sign == ' ' || sign == '-');

    double degs = double(std::get<1>(a_angle));   // Normally 0..359
    assert(0.0 <= degs && degs <= 359.0);

    double mins = double(std::get<2>(a_angle));   // Must be  0..59
    assert(0.0 <= mins && mins <= 59.0);

    double secs = std::get<3> (a_angle);          // Must be  0..60-
    assert(0.0 <= secs && secs <  60.0);

    // "Angle" is the same as "Angle_rad":
    Angle rads =
      To_Angle(Angle_deg   (double(degs))) +
      To_Angle(Angle_arcMin(double(mins))) +
      To_Angle(Angle_arcSec(secs));

    return (sign == '-') ? (-rads) : rads;
  }

  //=========================================================================//
  // "ToDMS":                                                                //
  //=========================================================================//
  // Conversion of an Angle in Radians into (sign, degs, arcMins, arcSecs), for
  // astronomical applications (eg Declination output):
  //
  constexpr inline std::tuple<double, Angle_deg, Angle_arcMin, Angle_arcSec>
    ToDMS(Angle a_angle)
  {
    // NB: Rounding is towards 0:
    double       sign = IsNeg(a_angle) ? -1.0 : IsPos(a_angle) ? 1.0 : 0.0;
    Angle_deg    degs = To_Angle_deg(Abs(a_angle));
    Angle_deg    d    = Floor(degs);
    Angle_arcMin mins = To_Angle_arcMin(degs - d);
    Angle_arcMin m    = Floor(mins);
    Angle_arcSec secs = To_Angle_arcSec(mins - m);

    assert(!(IsNeg(d) || IsNeg(m) || IsNeg(secs)));
    return std::make_tuple(sign, d, m, secs);
  }

  //=========================================================================//
  // "ToHMS":                                                                //
  //=========================================================================//
  // Conversion of an Angle in Radians into the range [0 .. 2*Pi), and then into
  // (hh, mm, ss), for astronomical applications (eg RightAscention output):
  //
  constexpr inline std::tuple<Angle_hh, Angle_mm, Angle_ss>
    ToHMS(Angle a_angle)
  {
    // Bring "a_angle" into [0 .. 2*Pi):
    Angle inRange = a_angle - Floor(a_angle / TwoPi<double>) * TwoPi<double>;

    // XXX: Is the following always true, even in the presence of rounding errs?
    assert(!IsNeg(inRange) && inRange < Angle(TwoPi<double>));

    Angle_hh  hours = To_Angle_hh(inRange);
    Angle_hh  hh    = Floor(hours);
    Angle_mm  mins  = To_Angle_mm(hours - hh);
    Angle_mm  mm    = Floor(mins);
    Angle_ss  ss    = To_Angle_ss(mins  - mm);

    return std::make_tuple(hh, mm, ss);
  }
}
// End namespace SpaceBallistics
