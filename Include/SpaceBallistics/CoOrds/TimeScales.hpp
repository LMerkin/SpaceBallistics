// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/CoOrds/TimeScales.hpp":								 //
//                     Non-ConstExpr TimeScales Functions										 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/PhysForces/DE440T.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "TDB" Ctor: From "TT":                                                  //
  //=========================================================================//
  TDB::TDB(TT a_tt)
    { *this = DE440T::TDBofTT(a_tt); }
}
// End namespace SpaceBallistics
