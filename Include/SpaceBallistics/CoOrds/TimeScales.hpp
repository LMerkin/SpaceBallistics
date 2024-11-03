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
  // "TT"  Ctor: From "TDB":                                                 //
  //=========================================================================//
  TT::TT(TDB const& a_tdb)
    { *this = DE440T::ToTT(a_tdb); }

  //=========================================================================//
  // "TDB" Ctor: From "TT":                                                  //
  //=========================================================================//
  TDB::TDB(TT const& a_tt)
    { *this = DE440T::ToTDB(a_tt); }
}
// End namespace SpaceBallistics
