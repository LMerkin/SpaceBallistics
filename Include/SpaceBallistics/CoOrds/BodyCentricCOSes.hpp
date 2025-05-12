// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/BodyCentricCOSes.hpp":              //
//        Translation of "Vector3D"s with Different Origins and Bodies       //
//===========================================================================//
// NB: These functions are NOT "constexpr", since they may use TT<->TDB conver-
// sions which are only available dynamically via DE440T:
//
#pragma  once
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/TimeScales.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // Translation: (BaryCEq, Earth) + (GeoCEqFix, Body) = (BaryCEq, Body):    //
  //=========================================================================//
  // XXX: This CANNOT be generalised to
  // (BaryEq, Body1) + (Body1Eq, Body2) = (BaryEq, Body2),
  // because BaryEq and Body1Eq COSes use DIFFERENT axes, unless Body1==Earth!
  //
  template<typename DQ,    Body B>
  Vector3D<DQ, BaryCEqCOS, B> operator+                     // Body  in BCRS
  (
    Vector3D<DQ, BaryCEqCOS,   Body::Earth> const& a_earth, // Earth in BCRS
    Vector3D<DQ, GeoCEqFixCOS, B>           const& a_body   // Body  in GCRS
  )
  {
    // XXX: Be careful with TimeStamps: "BaryCEqCOS" always uses TDB, whereas
    // "GeoCEqFix" uses TT. So unify the TimeStamps as TDB:
    TDB cosTSE(a_earth.GetCOSTS());  // Normally UnDef
    TDB cosTSB(a_body .GetCOSTS());  // TT->TDB using DE440T

    TDB vecTSE(a_earth.GetVecTS());
    TDB vecTSB(a_body .GetVecTS());  // TT->TDB using DE440T

    return
      Vector3D<DQ, BaryCEqCOS, B>
      { 
        UnifyTSs(cosTSE, cosTSB),
        UnifyTSs(vecTSE, vecTSB),
        a_earth.x() + a_body.x(),
        a_earth.y() + a_body.y(),
        a_earth.z() + a_body.z()
      };
  }

  //-------------------------------------------------------------------------//
  // Vector Addition is Commutative:                                         //
  //-------------------------------------------------------------------------//
  // So the following function must also be defined:
  //
  template<typename DQ,    Body B>
  Vector3D<DQ, BaryCEqCOS, B> operator+                     // Body  in BCRS
  (
    Vector3D<DQ, GeoCEqFixCOS, B>           const& a_body,  // Body  in GCRS
    Vector3D<DQ, BaryCEqCOS,   Body::Earth> const& a_earth  // Earth in BCRS
  )
  { return a_earth + a_body; }

  //=========================================================================//
  // Translation: (BaryCEq, Body) - (BaryCEq, Earth) = (GeoCEqFix, Body):    //
  //=========================================================================//
  // XXX: Again, this CANNOT be generalised to
  // (BaryCEq, Body1) - (BaryCEq, Body2) = (BodyCEq1, Body2):
  //
  template<typename DQ, Body B>
  Vector3D<DQ, GeoCEqFixCOS, B> operator-                   // Body  in GCRS
  (
    Vector3D<DQ, BaryCEqCOS, B>           const& a_body,    // Body  in BCRS
    Vector3D<DQ, BaryCEqCOS, Body::Earth> const& a_earth    // Earth in BCRS
  )
  {
    // NB: The resulting Vector uses TT, whereas the args use TDB:
    return
      Vector3D<DQ, GeoCEqFixCOS, B>
      {
        TT(UnifyTSs(a_earth.GetCOSTS(), a_body.GetCOSTS())),
        TT(UnifyTSs(a_earth.GetVecTS(), a_body.GetVecTS())),
        a_body.x() - a_earth.x(),
        a_body.y() - a_earth.y(),
        a_body.z() - a_earth.z()
      };
  }

  //-------------------------------------------------------------------------//
  // And other way round:                                                    //
  //-------------------------------------------------------------------------//
  template<typename DQ, Body B>
  Vector3D<DQ, GeoCEqFixCOS, B> operator-                   // - (Body in GCRS)
  (
    Vector3D<DQ, BaryCEqCOS, Body::Earth> const& a_earth,   // Earth in BCRS
    Vector3D<DQ, BaryCEqCOS, B>           const& a_body     // Body  in BCRS
  )
  { return - (a_body - a_earth); }

  //=========================================================================//
  // Translation: (BaryCEcl, Earth) + (GeoCEclFix, Body) = (BaryCEcl, Body): //
  //=========================================================================//
  // XXX: This CANNOT be generalised to
  // (BaryEcl, Body1) + (Body1Ecl, Body2) = (BaryEcl, Body2),
  // because BaryEcl and Body1Ecl COSes use DIFFERENT axes, unless Body1==Earth!
  //
  template<typename DQ,     Body B>
  Vector3D<DQ, BaryCEclCOS, B> operator+
  (
    Vector3D<DQ, BaryCEclCOS,   Body::Earth> const& a_earth,
    Vector3D<DQ, GeoCEclFixCOS, B>           const& a_body
  )
  {
    // XXX: Be careful with TimeStamps: "BaryCEclCOS" always uses TDB, whereas
    // "GeoCEclFix" uses TT. So unify the TimeStamps as TDB:
    TDB cosTSE(a_earth.GetCOSTS());  // Normally UnDef
    TDB cosTSB(a_body .GetCOSTS());  // TT->TDB using DE440T

    TDB vecTSE(a_earth.GetVecTS());
    TDB vecTSB(a_body .GetVecTS());  // TT->TDB using DE440T

    return
      Vector3D<DQ, BaryCEclCOS, B>
      { 
        UnifyTSs(cosTSE, cosTSB),
        UnifyTSs(vecTSE, vecTSB),
        a_earth.x() + a_body.x(),
        a_earth.y() + a_body.y(),
        a_earth.z() + a_body.z()
      };
  }

  //-------------------------------------------------------------------------//
  // Vector Addition is Commutative:                                         //
  //-------------------------------------------------------------------------//
  // So the following function must also be defined:
  //
  template<typename DQ,     Body B>
  Vector3D<DQ, BaryCEclCOS, B> operator+
  (
    Vector3D<DQ, GeoCEclFixCOS, B>           const& a_body,
    Vector3D<DQ, BaryCEclCOS,   Body::Earth> const& a_earth
  )
  { return a_earth + a_body; }

  //=========================================================================//
  // Translation: (BaryCEcl, Body) - (BaryCEcl, Earth) = (GeoCEclFix, Body): //
  //=========================================================================//
  // XXX: Again, this CANNOT be generalised to
  // (BaryCEcl, Body1) - (BaryCEcl, Body2) = (BodyCEcl1, Body2):
  //
  template<typename DQ,  Body B>
  Vector3D<DQ, GeoCEclFixCOS, B> operator-
  (
    Vector3D<DQ, BaryCEclCOS, B>           const& a_body,
    Vector3D<DQ, BaryCEclCOS, Body::Earth> const& a_earth
  )
  {
    // NB: The resulting Vector uses TT, whereas the args use TDB:
    return
      Vector3D<DQ, GeoCEclFixCOS, B>
      {
        TT(UnifyTSs(a_earth.GetCOSTS(), a_body.GetCOSTS())),
        TT(UnifyTSs(a_earth.GetVecTS(), a_body.GetVecTS())),
        a_body.x() - a_earth.x(),
        a_body.y() - a_earth.y(),
        a_body.z() - a_earth.z()
      };
  }

  //-------------------------------------------------------------------------//
  // And other way round:                                                    //
  //-------------------------------------------------------------------------//
  template<typename DQ,  Body B>
  Vector3D<DQ, GeoCEclFixCOS, B> operator-
  (
    Vector3D<DQ, BaryCEclCOS, Body::Earth> const& a_earth,
    Vector3D<DQ, BaryCEclCOS, B>           const& a_body
  )
  { return - (a_body - a_earth); }
}
// End namespace SpaceBallistics
