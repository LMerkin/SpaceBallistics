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
    // TODO: TimeInstant in some TimeScale...
    Mass                  m_fullMass;     // Curr full mass of a Stage
    Mass                  m_fuelMass;     // Curr mass of Fuel
    Mass                  m_oxidMass;     // Curr mass of Oxidiser
    MassRate              m_fullMassDot;  // Full Mass Time Derivative
    PosVEmb    <LVSCKind> m_com;          // Curr center of masses (x,  y,  z)
    VelVEmb    <LVSCKind> m_comDots;      // CoM  Embedded "Velocity"
    MoIVEmb    <LVSCKind> m_mois;         // Curr MoIs             (Jx, Jy, Jz)
    MoIRateVEmb<LVSCKind> m_moiDots;      // Curr MoIRates         (  similar )
    ForceVEmb  <LVSCKind> m_thrust;       // Curr engine thrust    (Fx, Fy, Fz)
    // TODO: Thrust Moment; AeroDynamic Forces and their Moments...
  };
}
// End namespace SpaceBallistics
