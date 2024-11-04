// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/EqEcl.hpp":										 //
//      Transformations between Ecliptical and Earth Equatorial Co-Ords      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"

#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // Vector Transformations between Equatorial and Ecliptical J2000.0 COSes: //
  //=========================================================================//
  // J2000.0 Obliquity of the Ecliptic to the ICRF Equator, from IAU76:
  //
  constexpr inline Angle     Obliq2000 = To_Angle(Angle_arcSec(84381.448));
  constexpr inline double CosObliq2000 = Cos(double(Obliq2000));
  constexpr inline double SinObliq2000 = Sin(double(Obliq2000));

  //-------------------------------------------------------------------------//
  // Eq -> Ecl:                                                              //
  //-------------------------------------------------------------------------//
  template<typename DQ,  typename EqCOS,  typename EclCOS>
  constexpr void  ToEcl
		(Vector3D<DQ, EqCOS> const& a_eq, Vector3D<DQ, EclCOS>* a_ecl)
  {
    // This transform makes sense for the following COSes only (as yet):
		static_assert
			((std::is_same_v<EqCOS,  BaryCentricEqCOS>   &&
				std::is_same_v<EclCOS, BaryCentricEclCOS>) ||
			 (std::is_same_v<EqCOS,  GeoCentricEqFixCOS> &&
				std::is_same_v<EclCOS, GeoCentricEclFixCOS>));
    assert(a_ecl != nullptr);

    a_ecl->x() =  a_eq.x();
    a_ecl->y() =  CosObliq2000 * a_eq.y() + SinObliq2000 * a_eq.z();
    a_ecl->z() = -SinObliq2000 * a_eq.y() + CosObliq2000 * a_eq.z();
  }
  
  template<typename DQ, typename EclCOS, typename EqCOS>
  constexpr Vector3D<DQ, EclCOS> ToEcl
           (Vector3D<DQ, EqCOS> const& a_eq)
  {
    Vector3D<DQ, BaryCentricEclCOS> ecl;
    ToEcl(a_eq, &ecl);
    return ecl;
  } 

  //-------------------------------------------------------------------------//
  // Ecl -> Eq:                                                              //
  //-------------------------------------------------------------------------//
  template<typename DQ,   typename EclCOS,  typename EqCOS>
  constexpr void  ToEq
    (Vector3D<DQ, EclCOS> const& a_ecl, Vector3D<DQ, EqCOS>* a_eq)
  {
    // This transform makes sense for the following COSes only (as yet):
		static_assert
			((std::is_same_v<EclCOS, BaryCentricEclCOS>   &&
				std::is_same_v<EqCOS,  BaryCentricEqCOS>) ||
			 (std::is_same_v<EclCOS, GeoCentricEclFixCOS> &&
				std::is_same_v<EqCOS,  GeoCentricEqFixCOS>));
    assert(a_eq != nullptr);

    a_eq->x() =  a_ecl.x();
    a_eq->y() =  CosObliq2000 * a_ecl.y() - SinObliq2000 * a_ecl.z();
    a_eq->z() =  SinObliq2000 * a_ecl.y() + CosObliq2000 * a_ecl.z();
  }

  template<typename DQ,  typename EclCOS, typename EqCOS>
  constexpr Vector3D<DQ, EclCOS>  ToEq
					 (Vector3D<DQ, EqCOS> const& a_ecl)
  {
    Vector3D<DQ, BaryCentricEclCOS> eq;
    ToEq(a_ecl, &eq);
    return eq;
  }
}
// End namespace SpaceBallistics
