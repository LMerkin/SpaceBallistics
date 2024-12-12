// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/PhysForces/DE440T-Data.h":                 //
//===========================================================================//
#pragma  once

namespace SpaceBallistics::DE440T::Bits
{
  //=========================================================================//
  // DE440T Data Dimensions:                                                 //
  //=========================================================================//
  // The Number of DE440T Data Records:
  constexpr inline int NR  = 5707;
  
  // Size of each record in "double"s:
  constexpr inline int ND  = 1122;

  // The Actual Data:
  extern double const Data[NR][ND];
}
// End namespace SpaceBallistics::DE440T::Bits
