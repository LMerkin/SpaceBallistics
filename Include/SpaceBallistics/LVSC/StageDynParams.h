// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/LVSC/StageDynParams.h":                    //
//               Masses, MoIs and Thrust of a Rocket Stage                   //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/LVSC/LVSC.h"
#include "SpaceBallistics/CoOrds/EmbeddedCOS.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "StageDynParams" Struct:                                                //
  //=========================================================================//
  // Masses, CoM, MoIs and the Thrust Vector (in the EMBEDDED Co-Ord System)
  // of an LV Stage at any given Time:
  //
  template<LVSC LVSCKind>
  struct StageDynParams
  {
    Mass                m_fullMass;  // Curr full mass of a Stage
    Mass                m_fuelMass;  // Curr mass of Fuel
    Mass                m_oxidMass;  // Curr mass if Oxidiser
    PosVEmb  <LVSCKind> m_com;       // Curr center of masses (x,  y,  z)
    MoITEmb  <LVSCKind> m_mois;      // Curr MoIs             (Jx, Jy, Jz)
    ForceVEmb<LVSCKind> m_thrust;    // Curr engine thrust    (Fx, Fy, Fz)
  };
}
// End namespace SpaceBallistics
