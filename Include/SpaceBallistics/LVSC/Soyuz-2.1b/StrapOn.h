// vim:ts=2:et
//===========================================================================//
//                "SpaceBallistics/LVSC/Soyuz-2.1b/StrapOn.h":               //
//  Mathematical Model of the "Soyuz-2.1b" StrapOn (any of "Blocks B,V,G,D") //
//===========================================================================//
#pragma  once 
#include "SpaceBallistics/ME/TrConeSpherSegm.hpp"
#include "SpaceBallistics/ME/ToricSegms.hpp"
#include "SpaceBallistics/LVSC/Propellants.h"
#include "SpaceBallistics/LVSC/StageDynParams.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Consts.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage2.h"
#include <cassert>

namespace SpaceBallistics
{
  namespace SC  = Soyuz21b_Consts;

  // All "MechElements" are instantiated with "LVSC::Soyuz21b":
  using     ME  = MechElement<LVSC::Soyuz21b>;
  using     PM  = PointMass  <LVSC::Soyuz21b>;
  using     TrC = TrCone     <LVSC::Soyuz21b>;
  using     SpS = SpherSegm  <LVSC::Soyuz21b>;
  using     Tor = ToricSegm  <LVSC::Soyuz21b>;
  using     S2  = Soyuz21b_Stage2;

  //=========================================================================//
  // "Soyuz21b_StrapOn" Class:                                               //
  //=========================================================================//
  template<unsigned I>
  class Soyuz21b_StrapOn
  {
  private:
    //=======================================================================//
    // No objects can be created of this class:                              //
    //=======================================================================//
    Soyuz21b_StrapOn() = delete;
    
  public:
    //=======================================================================//
    // Geometry:                                                             //
    //=======================================================================//
    static_assert(I <= 3);
    // I==0: Block B (+Y);
    // I==1: Block V (+Z);
    // I==2: Block G (-Y);
    // I==3: Block D (-Z):
    //
    // The Psi angle in the OYZ plane (see "RotationShell" for the definition):
    // Psi = Pi/2 * I:
    constexpr static double CosPsi = (I==0) ? 1.0 : (I==2) ? -1.0 : 0.0;
    constexpr static double SinPsi = (I==1) ? 1.0 : (I==3) ? -1.0 : 0.0;

    // The X-coord of the StrapOn top: Relative to MaxD of Stage2, ie,
    static_assert(S2::OxidTankUp.GetLow()[0] == S2::OxidTankLow.GetUp()[0]);
    constexpr static Len    TopX   = S2::OxidTankUp.GetLow()[0];
  };
}
// End namespace SpaceBallistivs
