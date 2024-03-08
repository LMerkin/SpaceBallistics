// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/StageDynParams.h":                    //
//               Masses, MoIs and Thrust of a Rocket Stage                   //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "StageDynParams" Struct:                                                //
  //=========================================================================//
  // Masses, CoM, MoIs and the Thrust Vector (in the EMBEDDED Co-Ord System)
  // of an LV Stage at any given Time:
  //
  struct StageDynParams
  {
    Mass  m_fullMass;   // Curr full mass of a Stage
    Mass  m_fuelMass;   // Curr mass of Fuel
    Mass  m_oxidMass;   // Curr mass if Oxidiser
    Len   m_com   [3];  // Curr center of masses (x,  y,  z)
    MoI   m_mois  [3];  // Curr MoIs             (Jx, Jy, Jz)
    Force m_thrust[3];  // Curr engine thrust    (Fx, Fy, Fz)
  };
}
