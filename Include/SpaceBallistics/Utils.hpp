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

  // Explicily allowing conversions of "Angle[_rad]" into "double":
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
    return Magnitude();
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
}
// End namespace SpaceBallistics
