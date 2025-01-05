// vim:ts=2:et
//===========================================================================//
//           "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.hpp":           //
//===========================================================================//
#pragma once
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b_Stage1_Booster::GetDynParams":                                //
  //=========================================================================//
  template<char Block>
  StageDynParams<LVSC::Soyuz21b>
  Soyuz21b_Stage1_Booster<Block>::GetDynParams
  (
    FlightTime       a_ft,
    Pressure         a_p,           // Curr Atmospheric Pressure
    ME::VelVE const& a_v,           // In the ECOS
    Angle_deg        a_verns_defl,  // Same angle for both Verns
    Angle_deg        a_aerofin_defl // AeroFin deflection angle
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    // We currently do not allow any times prior to LiftOff:
    assert(SC::LiftOffTime <= a_ft);

    // The atmospheric pressure is >= 0 obviously:
    assert(!IsNeg(a_p));

    //-----------------------------------------------------------------------//
    // Current Masses and Thrust:                                            //
    //-----------------------------------------------------------------------//

    StageDynParams<LVSC::Soyuz21b> res;
    return res;
  }
}
// End namespace SpaceBallistics
