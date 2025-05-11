// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/CoOrds/GeoCDynEqFixCOS.h":               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysEffects/EarthRotationModel.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/Maths/RotationMatrices.hpp"

namespace SpaceBallistics
{ 
  //=========================================================================//
  // "GeoCDynEqFixCOS":                                                      //
  //=========================================================================//
  // GeoCentric (NOT the general BodyCentric!) COS similar to GCRS, but the XY
  // plane and X axis are defined by the Dynamic True Equator  and the Dynamic
  // Ecliptic for the given Epoch, rather than J2000.0 as in GCRS (with a tiny
  // bias, actually).
  // NB: This type is not templated by the Epoch, because the latter may  be
  // Dynamic (with the Equator and the Equinox of a arbitrary date). Instead,
  // unlike other COS types, this type *does allow* construction of objects;
  // they contain the corresp EarthRotationMatrix:
  //
  class GeoCDynEqFixCOS
  {
  public:
    //-----------------------------------------------------------------------//
    // Consts and Types:                                                     //
    //-----------------------------------------------------------------------//
    constexpr static Body BaseBody       = Body::Earth;
    constexpr static bool HasFixedAxes   = true;
    constexpr static bool HasFixedOrigin = false;
    using  TimeScale = TT;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    // To convert vectors between this COS and GCRS, we need an Earth Axes Ori-
    // entation Matrix, so use a statically-computable "EarthRotationModel" for
    // that:
    EarthRotationModel const m_erm;

  public:
    //-----------------------------------------------------------------------//
    // Ctors and Accessors:                                                  //
    //-----------------------------------------------------------------------//
    // Default Ctor makes no sense:
    GeoCDynEqFixCOS() = delete;

    // Non-Default Ctor. NB: "m_erm" is constructed with Analytical Nutations
    // model (otherwise it could not be "constexpr"):
    //
    constexpr GeoCDynEqFixCOS(Time_jyr a_epoch)
    : m_erm(a_epoch)
    {}

    // Copy Ctor is auto-generated:
    constexpr GeoCDynEqFixCOS(GeoCDynEqFixCOS const&) = default;

    // Accessor:
    Time_jyr  GetEpoch() const { return m_erm.GetERMEpoch(); }

    //-----------------------------------------------------------------------//
    // "GeoCDynEqFixCOS" <-> GCRS Vector Conversions:                        //
    //-----------------------------------------------------------------------//
    // InvPN Matrix: Dyn  from GCRS
    // PN    Matrix: GCRS from Dyn
    //
    template<Body B, typename DQ>
    constexpr Vector3D<DQ, GeoCDynEqFixCOS, B> ToDynEquinox
      (Vector3D<DQ, GCRS, B> const& a_gcrs) const
    {
      Vector3D<DQ, GeoCDynEqFixCOS, B> res;
      m_erm.GetInvPN().MVMult(a_gcrs.GetArr(), res.GetArr());
      return res;
    }

    template<Body B, typename DQ>
    constexpr Vector3D<DQ, GCRS, B> ToGCRS
      (Vector3D<DQ, GeoCDynEqFixCOS, B> const& a_dyn) const
    {
      Vector3D<DQ, GCRS, B> res;
      m_erm.GetPN().MVMult(a_dyn.GetArr(), res.GetArr());
      return res;
    }
  };

  //-------------------------------------------------------------------------//
  // Vectors in "GeoCDynEqFixCOS":                                           //
  //-------------------------------------------------------------------------//
  template<Body B = Body::UNDEFINED>
  using DimLessV_GeoDynEqFix = DimLessV<GeoCDynEqFixCOS, B>;

  template<Body B = Body::UNDEFINED>
  using PosKV_GeoDynEqFix    = PosKV   <GeoCDynEqFixCOS, B>;

  template<Body B = Body::UNDEFINED>
  using VelKV_GeoDynEqFix    = VelKV   <GeoCDynEqFixCOS, B>;

}
// End namespace SpaceBallistics
