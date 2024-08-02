// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/LVSC/Soyuz-2.1b/LaunchPads.hpp":             //
//        Soyuz-2 (.1a, .1b) Pad Locations with Orientation Azimuths         //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/GeoLocations.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz2_LaunchPad" Class:                                               //
  //=========================================================================//
  // "Azimuth" is the azimuth (0..360 degs, from the North, Clock-Wise) of the
  // Main Axis of the Pad (ie pointing towards the Far Edge of the Pad).  Then
  // the initial Yaw (at the Launch Pad) can be taken as (Pi/4 - Azimuth):
  //
  class Soyuz2_LaunchPad: public Location<Body::Earth>
  {
    //-----------------------------------------------------------------------//
    // Data Fld(s):                                                          //
    //-----------------------------------------------------------------------//
  private:
    Angle   m_azimuth;
    Angle   m_yaw0;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // XXX: There is no much point in checking the range of the Azimuth:
    //
    constexpr Soyuz2_LaunchPad
    (
      Location<Body::Earth> const& a_location,
      Angle_deg                    a_azimuth  // Azimuth of the Main Pad Axis
    )
    : Location<Body::Earth>(a_location),
      m_azimuth            (To_Angle             (a_azimuth)),
      m_yaw0               (Angle(Pi_4<double>) - m_azimuth)
    {}

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Angle Azimuth() const { return m_azimuth; }
    constexpr Angle Yaw0   () const { return m_yaw0;    }
  };

  //=========================================================================//
  // Actual Pads:                                                            //
  //=========================================================================//
  constexpr static Soyuz2_LaunchPad Pad_Vostochny_1S  =
    Soyuz2_LaunchPad
    (
      Vostochny_1S,
      Location<Body::Earth>::GetAzimuth
        (128.33181_deg, 51.88161_deg, 128.33475_deg, 51.88435_deg)
    );

  constexpr static Soyuz2_LaunchPad Pad_Baykonur_31_6 =
    Soyuz2_LaunchPad
    (
      Baykonur_31_6,
      Location<Body::Earth>::GetAzimuth
        ( 63.5672_deg,  45.99428_deg,  63.56448_deg, 45.9959_deg)
    );

  constexpr static Soyuz2_LaunchPad Pad_Plesetsk_43_3 =
    Soyuz2_LaunchPad
    (
      Plesetsk_43_3,
      Location<Body::Earth>::GetAzimuth
        ( 40.45275_deg, 62.92641_deg,  40.45097_deg, 62.92689_deg)
    );

  constexpr static Soyuz2_LaunchPad Pad_Plesetsk_43_4 =
    Soyuz2_LaunchPad
    (
      Plesetsk_43_4,
      Location<Body::Earth>::GetAzimuth
        ( 40.4569_deg,  62.92725_deg,  40.45686_deg, 62.92803_deg)
    );
}
// End namespace SpaceBallistics
