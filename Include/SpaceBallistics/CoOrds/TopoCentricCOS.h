// vim:ts=2:et
//===========================================================================//
//                     "SpaceBallistics/CoOrds/CoOrds.h":                    //
//                   Co-Ordinate Systems and State Vectors                   //
//===========================================================================//
#pragma once
#include "SpaceBallistics/CoOrds/Locations.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // Co-Ord Systems:                                                         //
  //=========================================================================//
  //=========================================================================//
  // "TopoCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: A point on the Earth surface given by "L"
  // Axes  : (X=East, Y=North, Z=Zenith)
  // NB    : This type  only stands for itself;
  //         no objects of it are to be created:
  //
  template<Location_WGS84 const* L>
  class TopoCentricCOS
  {
    TopoCentricCOS() = delete;
  };
}
// End namespave SpaceBallistics
