// vim:ts=2:et
//===========================================================================//
//                               "Soyuz21b.h":                               //
//          Mathematical Model of the "Soyuz-2.1b" Launch Vehicle            //
//===========================================================================//
// Includes Blocks A, B-D, I, InterStage and Fairing, but does NOT include the
// 4th stage (eg "Fregate"), the payload adapter / dispenser  and  the Payload
// itself:
//
#pragma  once
#include "SpaceBallistics/Types.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Soyuz21b" Class:                                                       //
  //=========================================================================//
  // For MoI computations, the following local co-ord system is used for con-
  // venience (so that most X-coords are positive):
  // (*) The OX axis is the main axis of the rocket. The OX positive direction
  //     is towards the TAIL (REAL of the LV). The origin O is the upper base
  //     of the InterStage (junction plane with the Fairing). That is, X < 0
  //     for the optional Stage4, payload adapter/dispenser,  payload itself
  //     and the fairing, and X >= 0 for the Integstage and Stages 3, 2, 1.
  // (*) The OY and OZ axes are such that the OXY and OXZ planes pass through
  //     the symmetry axes of the corresp opposite strap-on boosters (Blocks
  //     B, V, G, D -- Stage 1), and OXYZ is a right-oriented co-ords system.
  // XXX:
  // (*) All computations here are *APPROXIMATE*, using the thin shell model
  //     and point-mass (for the engines) models. The results are normalised
  //     to the empty masses of the rocket stages which are known reasonably
  //     accurately.
  //
  class Soyuz21b
  {
  public:
    //=======================================================================//
    // Static Data:                                                          //
    //=======================================================================//
    // XXX: The following consts would better be made "constexpr"s, but we cur-
    // rently cannot do that because their computation involves non-"constexpr"
    // funcs. It would also be possible to make them "static",  but then their
    // initialisation becomes more difficult. So they are made member flds:
    //
    Mass const  m_stage3EmptyMass;      // With InterStage and Aft
    MoI  const  m_stage3EmptyMoIY;
    Mass const  m_stage3NoAftEmptyMass; // Same but w/o Aft
    MoI  const  m_stage3NoAftEmptyMoIY;

    //=======================================================================//
    // Default and Copy Ctors:                                               //
    //=======================================================================//
    Soyuz21b();
    Soyuz21b(Soyuz21b const& a_right);
  };
}
