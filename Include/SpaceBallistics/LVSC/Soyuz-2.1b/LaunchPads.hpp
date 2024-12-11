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
  // "Azimuth" is the Pad Orientation Azimuth, is the angle (0..360 degs, from
  // the North, Clock-Wise) of the Main Axis of the Pad (ie pointing towards
  // the Far Edge of the Pad).
  // Then the initial LV Yaw  (at  the Launch Pad, counter-clock-wise)  can be
  // taken as (Pi/4 - Azimuth), in the "TopoCRotCOS".
  // However, in the "TopoCLaunchCOS" the initial Yaw depends on the Launch Azi-
  // muth, not on the Pad Orientation Azimuth...
  //
  class Soyuz2_LaunchPad: public Location<Body::Earth>
  {
    //-----------------------------------------------------------------------//
    // Data Fld(s):                                                          //
    //-----------------------------------------------------------------------//
  private:
    Angle   m_padAzimuth;

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
      m_padAzimuth         (To_Angle(a_azimuth))
    {}

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    constexpr Angle PadAzimuth() const { return m_padAzimuth; }
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
