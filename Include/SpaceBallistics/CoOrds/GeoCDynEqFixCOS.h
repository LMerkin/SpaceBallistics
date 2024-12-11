// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/CoOrds/GeoCDynEqFixCOS.h":               //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysForces/EarthRotationModel.h"

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
  // XXX: "EpochY" is just an "int", not "Time_jyr", because the latter is not a
  // structural type (contains a private member), and cannot therefore be used
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
    // entation Matrix, so use a statocally-computable "EarthRottaionModel" for
    // that:
    constexpr static EarthRotationModel ERM = EarthRotationModel(EpochY);
  };
}
// End namespace SpaceBallistics
