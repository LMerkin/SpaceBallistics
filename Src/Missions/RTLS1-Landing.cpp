// vim:ts=2:et
//===========================================================================//
//                     "Src/Missions/RTLS1-Landing.cpp":                     //
//===========================================================================//
#include "SpaceBallistics/Missions/RTLS1.h"

namespace SpaceBallistics
{
//===========================================================================//
// "SetLandBurnParams":                                                      //
//===========================================================================//
// This function is constructed from the data obtained using
// "CalibrateLandBurnParams".
// It maps the observed velocities (just after crossing the "MaxLandBurnH"
// boundary) to the optimal "LandBurnH" and "LandBurnGamma"  which  should
// provide the soft landing (with Velocity and Acceleration constraints sa-
// tisfied):
//
//
void RTLS1::SetLandBurnParams
(
  LenK UNUSED_PARAM(a_ref_h),
  VelK UNUSED_PARAM(a_ref_Vr),
  VelK UNUSED_PARAM(a_ref_Vhor)
)
{
  // XXX: For the moment, it is a placeholder only: No LandBurn actually occurs,
  // we just hit the Earth surface:
  m_landBurnH     = 0.0_km;
  m_landBurnGamma = 0.0;
}
}
