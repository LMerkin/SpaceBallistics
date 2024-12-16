// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/CoOrds/GeoCDynEqFixCOS.h":               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysForces/EarthRotationModel.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"

namespace SpaceBallistics
{ 
  //=========================================================================//
  // "GeoCDynEqFixCOS":                                                      //
  //=========================================================================//
  // GeoCentric (NOT the general BodyCentric!) COS similar to GCRS,  but the XY
  // plane and X axis are defined  by the Dynamical True Equator  and  the Mean
  // Ecliptics for the given Epoch Year, rather than J2000.0 (with a tiny bias)
  // as in GCRS.
  // XXX: Mostly intended for testing purposes; production use in inconvenient
  // because tis type is templated by the Epoch (formally in TT), so that mul-
  // tiple Epochs would need to be created statically.
  // XXX: "EpochY" is just an "int", not "Time_jyr", because the latter is not
  // a structural type (contains a private member), and cannot therefore be used
  // as a template parameter:
  //
  template<int EpochY>
  struct GeoCDynEqFixCOS
  {
    constexpr static Body BaseBody       = Body::Earth;
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = false;
    using  TimeScale = TT;

    // This struct stands for itself; no objects of it can be created:
    GeoCDynEqFixCOS() = delete;

    // To convert vectors between this COS and GCRS, we need an Earth Axes Ori-
    // entation Matrix, so use a statically-computable "EarthRotationModel" for
    // that:
    constexpr static EarthRotationModel ERM = EarthRotationModel(EpochY);
  };

  //-------------------------------------------------------------------------//
  // Vectors in "GeoCDynEqFixCOS":                                           //
  //-------------------------------------------------------------------------//
  template<int EpochY, Body B = Body::UNDEFINED>
  using PosKV_GeoDynEqFix = PosKV<GeoCDynEqFixCOS<EpochY>, B>;

  template<int EpochY, Body B = Body::UNDEFINED>
  using VelKV_GeoDynEqFix = VelKV<GeoCDynEqFixCOS<EpochY>, B>;

  //-------------------------------------------------------------------------//
  // "GeoCDynEqFixCOS" <-> GCRS Vector Conversions:                          //
  //-------------------------------------------------------------------------//
  // InvPN Matrix: Dyn  from GCRS
  // PN    Matrix: GCRS from Dyn
  //
  template<int EpochY,  Body B, typename DQ>
  constexpr Vector3D<DQ, GeoCDynEqFixCOS<EpochY>, B> ToDynEquinox
    (Vector3D<DQ, GCRS, B> const& a_gcrs)
  {
    Vector3D<DQ, GeoCDynEqFixCOS<EpochY>, B> res;
    GeoCDynEqFixCOS<EpochY>::ERM.GetInvPN().MVMult
      (a_gcrs.GetArr(), res.GetArr());
    return res;
  }

  template<int EpochY, Body B, typename DQ>
  constexpr Vector3D<DQ, GCRS, B> ToGCRS
    (Vector3D<DQ, GeoCDynEqFixCOS<EpochY>, B> const& a_dyn)
  {
    Vector3D<DQ, GCRS, B> res;
    GeoCDynEqFixCOS<EpochY>::ERM.GetPN().MVMult
      (a_dyn.GetArr(), res.GetArr());
    return res;
  }
}
// End namespace SpaceBallistics
