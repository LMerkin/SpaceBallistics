// vim:ts=2:et
//===========================================================================//
//                         "Soyuz2_LaunchPads.hpp":                          //
//                         Locations with Azimuths                           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz2_LaunchPad" Class:                                               //
  //=========================================================================//
  // "Azimuth" is the azimuth (0..360 degs, from the North, Clock-Wise) of the
  // Main Axis of the Pad (ie pointing towards the Far Edge of the Pad).  Then
  // the initial Yaw (at the Launch Pad) can be taken as (Pi/4 - Azimuth):
  //
  class Soyuz2_LaunchPad: public Location_WGS84
  {
    //-----------------------------------------------------------------------//
    // Data Fld(s):                                                          //
    //-----------------------------------------------------------------------//
  private:
    Angle   m_azimuth;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // XXX: There is no much point in checking the range of the Azimuth:
    //
    constexpr Soyuz2_LaunchPad
    (
      Location_WGS84 a_location,
      Angle_deg      a_azimuth   // Azimuth of the Main Pad Axis
    )
    : Location_WGS84(a_location),
      m_azimuth     (To_Angle_rad(a_azimuth))
    {}

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Angle Azimuth() const { return m_azimuth; }
    constexpr Angle Yaw0   () const { return Pi_4<double> - m_azimuth; }
  };

  //=========================================================================//
  // Actual Pads:                                                            //
  //=========================================================================//
  constexpr static Soyuz2_LaunchPad Pad_Vostochny_1S  =
    Soyuz2_LaunchPad(Vostochny_1S, );

  constexpr static Soyuz2_LaunchPad Pad_Baykonur_31_6 =
    Soyuz2_LaunchPad(Baykonur_31_6, );

  constexpr static Soyuz2_LaunchPad Pad_Plesetsk_43_3 =
    Soyuz2_LaunchPad(Plesetsk_43_3, );

  constexpr static Soyuz2_LaunchPad Pad_Plesetsk_43_4 =
    Soyuz2_LaunchPad(Plesetsk_43_4, )
}
