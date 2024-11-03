// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/BaryCentricCOSes.h":                 //
//                    The Inertial Co-Ordinate System                        //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "BaryCentricEclCOS" Class:                                              //
  //=========================================================================//
  // BaryCentric Ecliptical COS:
  // Origin: Solar System BaryCenter
  // Axes  :
  //   X       : ICRF Equinox of J2000.0
  //   XY plane: Ecliptic     of J2000.0
  //             (ie NOT ICRF, which is Equatorial, and NOT the Laplace plane
  //             of the Solar System, which is the plane orthogonal to the
  //             MoI vector
  //   Y, Z    : Derived from the above
  // Obviously, using TDB as the assiciated TimeScale.
  // This definition is consistent with the conventions of the JPL Horizon eph-
  // emerides:
  //
  class  TDB;

  struct BaryCentricEclCOS
  {
    using TimeScale  = TDB;

    BaryCentricEclCOS() = delete;   // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  using PosKVBEcl = PosKV<BaryCentricEclCOS>;
  using VelKVBEcl = VelKV<BaryCentricEclCOS>;

  // XXX: Currently no need to consider other Vectors in this COS yet...

  //=========================================================================//
  // "BaryCentricEqCOS" Class:                                               //
  //=========================================================================//
  // BaryCentric Equatorial COS:
  // As "BaryCentricEclCOS" above, but the XY plane is the Earth Equator of
  // J2000.0. Thus, the axes of this COS coincide with the ICRF. Might   be
  // useful eg for integration of GeoCentric trajectories, with perturbations
  // from the Moon and the Sun:
  //
  struct BaryCentricEqCOS
  {
    using TimeScale    = TDB;
    BaryCentricEqCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  using PosKVBEq  = PosKV<BaryCentricEqCOS>;
  using VelKVBEq  = VelKV<BaryCentricEqCOS>;

  // XXX: Currently no need to consider other Vectors in this COS yet...

  //=========================================================================//
  // Transformations of Vectors between Equatorial and Ecliptic Co-Ords:     //
  //=========================================================================//
  // J2000.0 Obliquity of the Ecliptic to the ICRF Equator, from IAU76:
  //
  constexpr inline Angle     Obliq2000 = To_Angle(Angle_arcSec(84381.448));
  constexpr inline double CosObliq2000 = Cos(double(Obliq2000));
  constexpr inline double SinObliq2000 = Sin(double(Obliq2000));

  //-------------------------------------------------------------------------//
  // Eq -> Ecl:                                                              //
  //-------------------------------------------------------------------------//
  template<typename DQ>
  constexpr void  ToEcl
    (Vector3D<DQ, BaryCentricEqCOS> const& a_eq,
     Vector3D<DQ, BaryCentricEclCOS>*      a_ecl)
  {
    assert(a_ecl != nullptr);
    a_ecl->x() =  a_eq.x();
    a_ecl->y() =  CosObliq2000 * a_eq.y() + SinObliq2000 * a_eq.z();
    a_ecl->z() = -SinObliq2000 * a_eq.y() + CosObliq2000 * a_eq.z();
  }

  template<typename DQ>
  constexpr Vector3D<DQ, BaryCentricEclCOS> ToEcl
           (Vector3D<DQ, BaryCentricEqCOS> const& a_eq)
  {
    Vector3D<DQ, BaryCentricEclCOS> ecl;
    ToEcl(a_eq,  &ecl);
    return ecl;
  }

  //-------------------------------------------------------------------------//
  // Ecl -> Eq:                                                              //
  //-------------------------------------------------------------------------//
  template<typename DQ>
  constexpr void  ToEq
    (Vector3D<DQ, BaryCentricEclCOS> const& a_ecl,
     Vector3D<DQ, BaryCentricEqCOS>*        a_eq)
  {
    assert(a_eq != nullptr);
    a_eq->x() =  a_ecl.x();
    a_eq->y() =  CosObliq2000 * a_ecl.y() - SinObliq2000 * a_ecl.z();
    a_eq->z() =  SinObliq2000 * a_ecl.y() + CosObliq2000 * a_ecl.z();
  }

  template<typename DQ>
  constexpr Vector3D<DQ, BaryCentricEclCOS> ToEq
           (Vector3D<DQ, BaryCentricEqCOS> const& a_ecl)
  {
    Vector3D<DQ, BaryCentricEclCOS> eq;
    ToEq(a_ecl,  &eq);
    return eq;
  }
}
// End namespace SpaceBallistics
