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
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // Vector Transformations between Equatorial and Ecliptical ~J2000 COSes:  //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // For the Earth Ecliptic, Equator and Equinox of ~J2000.0:                //
  //-------------------------------------------------------------------------//
  // J2000.0 Obliquity of Ecliptic to the ICRS Equator, from the IAU76 model:
  //
  constexpr inline Angle     Obliq2000 = To_Angle(EarthRotationModel::Eps0);
  constexpr inline double CosObliq2000 = Cos(double(Obliq2000));
  constexpr inline double SinObliq2000 = Sin(double(Obliq2000));

  //-------------------------------------------------------------------------//
  // BaryCentrocEqCOS   -> BaryCentricEclCOS:                                //
  // GeoCentricEqFixCOS -> GeoCentricEclFixCOS:                              //
  //-------------------------------------------------------------------------//
  template
  <typename DQ, typename EqCOS, typename EclCOS, Body B = Body::UNDEFINED>
  constexpr void  ToEcl
  (
    Vector3D<DQ, EqCOS,  B> const& a_eq,
    Vector3D<DQ, EclCOS, B>*       a_ecl
  )
  requires((std::is_same_v<EqCOS,  BaryCentricEqCOS>   &&
            std::is_same_v<EclCOS, BaryCentricEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCentricEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCentricEclFixCOS>))
  {
    assert(a_ecl != nullptr);
    a_ecl->x() =  a_eq.x();
    a_ecl->y() =  CosObliq2000 * a_eq.y() + SinObliq2000 * a_eq.z();
    a_ecl->z() = -SinObliq2000 * a_eq.y() + CosObliq2000 * a_eq.z();
  }

  template
  <typename  DQ, typename EqCOS, typename EclCOS, Body B = Body::UNDEFINED>
  constexpr Vector3D<DQ, EclCOS, B> ToEcl(Vector3D<DQ, EqCOS, B> const& a_eq)
  requires((std::is_same_v<EqCOS,  BaryCentricEqCOS>   &&
            std::is_same_v<EclCOS, BaryCentricEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCentricEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCentricEclFixCOS>))
  {
    Vector3D<DQ, EclCOS, B> ecl;
    ToEcl(a_eq, &ecl);
    return       ecl;
  }

  //-------------------------------------------------------------------------//
  // BaryCentricEclCOS   -> BaryCentricEqCOS:                                //
  // GeoCentricEclFixCOS -> GeoCentricEclFixCOS:                             //
  //-------------------------------------------------------------------------//
  template
  <typename DQ, typename EclCOS, typename EqCOS, Body B = Body::UNDEFINED>
  constexpr void ToEq
  (
    Vector3D<DQ, EclCOS, B> const& a_ecl,
    Vector3D<DQ, EqCOS,  B>*       a_eq
  )
  requires((std::is_same_v<EqCOS,  BaryCentricEqCOS>   &&
            std::is_same_v<EclCOS, BaryCentricEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCentricEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCentricEclFixCOS>))
  {
    assert(a_eq != nullptr);
    a_eq->x() =  a_ecl.x();
    a_eq->y() =  CosObliq2000 * a_ecl.y() - SinObliq2000 * a_ecl.z();
    a_eq->z() =  SinObliq2000 * a_ecl.y() + CosObliq2000 * a_ecl.z();
  }

  template
  <typename  DQ, typename  EclCOS, typename EqCOS, Body B = Body::UNDEFINED>
  constexpr Vector3D<DQ, EqCOS, B> ToEq(Vector3D<DQ, EclCOS, B> const& a_ecl)
  requires((std::is_same_v<EqCOS,  BaryCentricEqCOS>   &&
            std::is_same_v<EclCOS, BaryCentricEclCOS>) ||
           (std::is_same_v<EqCOS,  GeoCentricEqFixCOS> &&
            std::is_same_v<EclCOS, GeoCentricEclFixCOS>))
  {
    Vector3D<DQ, EqCOS, B> eq;
    ToEq(a_ecl, &eq);
    return       eq;
  }
}
// End namespace SpaceBallistics
