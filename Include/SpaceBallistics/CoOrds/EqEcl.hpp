// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/EqEcl.hpp":                     //
//      Transformations between Ecliptical and Earth Equatorial Co-Ords      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // Vector Transformations between Equatorial and Ecliptical ~J2000 COSes:  //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // For the Earth Ecliptic, Equator and Equinox of ~J2000.0:                //
  //-------------------------------------------------------------------------//
  // XXX: In the "EarthRotationModel" and here we assume that the axes of Eq2000
  // are the same as those of ICRS/BCRS/GCRS, ie the bias is disregarded:
  //
  // J2000.0 Obliquity of Ecliptic to the ICRS Equator:
  constexpr inline Angle  EObliqJ2000 = To_Angle(EarthRotationModel::Eps0);

  // Rotation Matrices:
  constexpr inline Mtx33  ToEclJ2000  = Mtx33::MkR1(EObliqJ2000);
  constexpr inline Mtx33  ToEqJ2000   = ToEclJ2000.Transpose();

  //-------------------------------------------------------------------------//
  // BaryCEqCOS   -> BaryCEclCOS:                                            //
  // GeoCEqFixCOS -> GeoCEclFixCOS:                                          //
  //-------------------------------------------------------------------------//
  // "EclCOS" typically needs to be specified explicitly, other template params
  // may be inferred:
  template
  <typename EclCOS, typename DQ, typename EqCOS, Body B = Body::UNDEFINED>
  constexpr Vector3D<DQ, EclCOS, B> ToEcl(Vector3D<DQ, EqCOS, B> const& a_eq)
  requires((std::is_same_v<EqCOS,  BaryCEqCOS>   &&
            std::is_same_v<EclCOS, BaryCEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCEclFixCOS>))
  {
    Vector3D<DQ, EclCOS, B> ecl;
    ToEclJ2000.MVMult(a_eq.GetArr(), ecl.GetArr());
    return ecl;
  }

  //-------------------------------------------------------------------------//
  // BaryCEclCOS   -> BaryCEqCOS:                                            //
  // GeoCEclFixCOS -> GeoCEclFixCOS:                                         //
  //-------------------------------------------------------------------------//
  // "EqCOS" typically needs to be specified explicitly, other template params
  // may be inferred:
  template
  <typename EqCOS, typename  DQ, typename  EclCOS, Body B = Body::UNDEFINED>
  constexpr Vector3D<DQ, EqCOS, B> ToEq(Vector3D<DQ, EclCOS, B> const& a_ecl)
  requires((std::is_same_v<EqCOS,  BaryCEqCOS>   &&
            std::is_same_v<EclCOS, BaryCEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCEclFixCOS>))
  {
    Vector3D<DQ, EqCOS, B> eq;
    ToEqJ2000.MVMult(a_ecl.GetArr(), eq.GetArr());
    return eq;
  }
}
// End namespace SpaceBallistics
