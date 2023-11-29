// vim:ts=2:et
//===========================================================================//
//                                 "Soyuz21b.cpp":                           //
//          Mathematical Model of the "Soyuz-2.1b" Launch Vehicle            //
//===========================================================================//
#include "SpaceBallistics/Soyuz21b.h"
#include "SpaceBallistics/Soyuz21b.hpp"
#include "SpaceBallistics/ConstrElement.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // Default Ctor:                                                           //
  //=========================================================================//
  Soyuz21b::Soyuz21b()
  {
    using namespace Soyuz21b_Consts;

    //-----------------------------------------------------------------------//
    // IntegStage and Stage 3 ("Block I"):                                   //
    //-----------------------------------------------------------------------//
    // Compute the empty Mass, MoI and CoM with and w/o the Aft:
    //
    // "shells3" is Stage3 incl the Aft but w/o the Engine (which is modeled as
    // a PointMass);
    // Begin with the over-all cylindrical shell of Stage3:
    ConstrElement st3 =
      CE::ShellTrCone(Stage3X0, 0.0_m,   0.0_m,   0.0,
                      Stage3D,  Stage3D, Stage3H, 0.0_kg);

    // Stage3FuelTank: UpperSpherSegm (neg-facing),
    //                 LowerSpherSegm (pos-facing):
    // NB: The cylindrical (middle) section of Stage3FuelTank is NOT included
    // in the calculations as it is accounted for in the "bigCyl3":
    st3 +=
      CE::ShellSpherSegm
        (false,   Stage3FuelTankUpX0, 0.0_m, 0.0_m, 0.0,
         Stage3D, Stage3FuelTankUpH,  0.0_kg);
    st3 +=
      CE::ShellSpherSegm
        (true,    Stage3FuelTankLoX0, 0.0_m, 0.0_m, 0.0,
         Stage3D, Stage3FuelTankLoH,  0.0_kg);

    // Stage3 Oxidiser Tank: Similar to the Fuel Tank:
    // The cylindrical (middle) section is NOT included for the above reason:
    st3 +=
      CE::ShellSpherSegm
        (false,   Stage3OxidTankUpX0, 0.0_m, 0.0_m, 0.0,
         Stage3D, Stage3OxidTankUpH,  0.0_kg);
    st3 +=
      CE::ShellSpherSegm
        (true,    Stage3OxidTankLoX0, 0.0_m, 0.0_m, 0.0,
         Stage3D, Stage3OxidTankLoH,  0.0_kg);

    // Now set the actual mass of "st3" (which corresponds to Stage3 shells w/o
    // the Engine), and pro-rate the MoI:
    st3.SetMass(Stage3EmptyMass - Stage3EngMass);

    // Only now we can add the InterStage with its own mass:
    st3 +=
      CE::ShellTrCone(InterStageX0, 0.0_m, 0.0_m, 0.0,
                      InterStageD0, InterStageD1, InterStageH,
                      InterStageMass);

    // The jettisonable Aft section on its own (with initially unknown mass):
    ConstrElement aft3 =
      CE::ShellTrCone(Stage3OxidTankLoX0, 0.0_m, 0.0_m, 0.0,
                      Stage3D,   Stage3D, Stage3AftH,   0.0_kg);

    // Now, assuming the surface density of the Aft section is the same as for
    // "st3" over-all, set its mass:
    Mass aft3Mass =
      double(aft3.GetSurfArea() / st3.GetSurfArea()) * st3.GetMass();
    aft3.SetMass(aft3Mass);

    // Then add the Stage3 Engine as a Point Mass:
    st3 += CE::PointMass(Stage3EngX0, 0.0_m, 0.0_m, Stage3EngMass);

    // Set the computed Stage3 Consts:
    const_cast<Mass&>(m_stage3EmptyMass)      = st3.GetMass();
    const_cast<MoI&> (m_stage3EmptyMoIY)      = st3.GetMoIs()[1];

    const_cast<Mass&>(m_stage3NoAftEmptyMass) =
      m_stage3EmptyMass - aft3Mass;
    assert(IsPos(m_stage3NoAftEmptyMass));

    const_cast<MoI&> (m_stage3NoAftEmptyMoIY) =
      m_stage3EmptyMoIY - aft3.GetMoIs()[1];
    assert(IsPos(m_stage3NoAftEmptyMoIY));
  }
}
