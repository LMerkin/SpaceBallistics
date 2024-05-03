// vim:ts=2:et
//===========================================================================//
//                        "SpaceBallistics/LVSC/LVSC.h":                     //
//   Defines a Nomenclature of all Launch Vehicles and SpaceCraft Modeled    //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "LVSC" Enum Class:                                                      //
  //=========================================================================//
  enum class LVSC: int
  {
    BoilerPlate0,  // Dummy PayLoad
    Soyuz21b       // Soyuz-2.1b LV (w/o the PayLoad)
    // Other LVSC Kinds to come...
  };
}
// End namespace SpaceBallistics
