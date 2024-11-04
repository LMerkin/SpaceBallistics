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
  //             angular momentum vector of the whole Solar System)
  //   Y, Z    : Derived from the above
  // Obviously, using TDB as the associated TimeScale.
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
}
// End namespace SpaceBallistics
