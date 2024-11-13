// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/CoOrds/BaryCentricCOSes.h":                 //
//                    The Inertial Co-Ordinate System                        //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  // Fwd Decls:
  class  TDB;

  //=========================================================================//
  // "BaryCentricEclCOS" Struct:                                             //
  //=========================================================================//
  // BaryCentric Ecliptical COS:
  // Origin: Solar System BaryCenter
  // Axes  :
  //   X       : ICRS Equinox of J2000.0
  //   XY plane: Ecliptic     of J2000.0
  //             (ie NOT ICRS, which is Equatorial, and NOT the Laplace plane
  //             of the Solar System, which is the plane orthogonal to the
  //             angular momentum vector of the whole Solar System)
  //   Y, Z    : Derived from the above
  // Obviously, using TDB as the associated TimeScale.
  // This definition is consistent with the conventions of the JPL Horizon eph-
  // emerides:
  //
  struct BaryCentricEclCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCentricEclCOS() = delete;  // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  using PosKVBarEcl = PosKV<BaryCentricEclCOS>;
  using VelKVBarEcl = VelKV<BaryCentricEclCOS>;

  // XXX: Currently no need to consider other Vectors in this COS yet...

  //=========================================================================//
  // "BaryCentricEqCOS" Struct:                                              //
  //=========================================================================//
  // BaryCentric Equatorial COS:
  // As "BaryCentricEclCOS" above, but the XY plane is the Earth Equator of
  // J2000.0. Thus, the axes of this COS coincide with the ICRS. Might   be
  // useful eg for integration of GeoCentric trajectories, with perturbations
  // from the Moon and the Sun:
  //
  struct BaryCentricEqCOS
  {
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = true;
    using TimeScale                      = TDB;
    BaryCentricEqCOS() = delete;   // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors and Tensors in this COS:                  //
  //-------------------------------------------------------------------------//
  using PosKVBarEq  = PosKV<BaryCentricEqCOS>;
  using VelKVBarEq  = VelKV<BaryCentricEqCOS>;

  // XXX: Currently no need to consider other Vectors in this COS yet...
}
// End namespace SpaceBallistics
