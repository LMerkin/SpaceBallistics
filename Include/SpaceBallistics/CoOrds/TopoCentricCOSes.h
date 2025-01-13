// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/TopoCentrticCOSes.h":               //
//                      Body-TopoCentric Co-Ord Systems                      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/Locations.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "TopoCEqFixCOS" Struct:                                                 //
  //=========================================================================//
  // The Axes are as in "BodyCEqFixCOS" (in particular, for Earth they are ICRS/
  // GCRS), but the Origin is a Location ("L") on the Body, NOT the Body Center.
  // Suitable for computation of TopoCentric Ephemerides in astronomical appls.
  // The Origin of this COS is rotating (and otherwise moving) in the inertial
  // space, but the Axes are fixed.
  // IMPORTANT: Similar to BodyCentric COSes, TopoCentric COSes are SnapShots,
  // so eg Location's Velocity and Acceleration  are NOT  identical 0s in the
  // resp TopoCentric COS. The TimeStamp of such a SnapShot cannot be installed
  // in the COS type itself (because it is typically known at run-time only);
  // rather, it is inslalled in the corresp "Vector3D"s:
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCEqFixCOS
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = true;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCEqFixCOS() = delete; // No objs of this struct are to be created
  };

  // In particular, for Earth we call it "TopCGCRS":
  template<Location<Body::Earth> const* L>
  using TopoC_GCRS  = TopoCEqFixCOS<Body::Earth, L>;

  //-------------------------------------------------------------------------//
  // Position and Velocity Vectors in this COS:                              //
  //-------------------------------------------------------------------------//
  // We probably do not need any other DQs here:
  //
  template<Body BBody,  Location<BBody> const* L, Body B = Body::UNDEFINED>
  using DimLessVTopEqFix = DimLessV<TopoCEqFixCOS<BBody, L>, B>;

  template<Body BBody,  Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopEqFix    = PosKV   <TopoCEqFixCOS<BBody, L>, B>;

  template<Body BBody,  Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopEqFix    = VelKV   <TopoCEqFixCOS<BBody, L>, B>;

  // Aliases for Vectors in the "TopoC_GCRS":
  template<Location<Body::Earth> const* L,        Body B = Body::UNDEFINED>
  using DimLessTopGCRS   = DimLessV<TopoC_GCRS<L>,  B>;

  template<Location<Body::Earth> const* L,        Body B = Body::UNDEFINED>
  using PosKVTopGCRS     = PosKV   <TopoC_GCRS<L>,  B>;

  template<Location<Body::Earth> const* L,        Body B = Body::UNDEFINED>
  using VelKVTopGCRS     = VelKV   <TopoC_GCRS<L>,  B>;

  //=========================================================================//
  // "TopoCRotCOS_ENZ" Struct:                                               //
  //=========================================================================//
  // Origin   : A point on the Body Surface given by "L";
  // Axes     : X=LocalEast   (tangential to the Parallel);
  //            Y=LocalNorth  (tangential to the Meridian);
  //            Z=LocalZenith (normal to the Ellipsoid Surface);
  // TimeScale: TT for Earth, TDB otherwise;
  // This COS is HIGHLY NON-INERTIAL: its Origin is rotating (and otherwise mov-
  // ing) in the inertial space, and the axes are rotating as well.
  // Suitable eg for modeling of Launch Ops and for Ground Tracks!
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCRotCOS_ENZ
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCRotCOS_ENZ() = delete;
    // No objs of this struct are to be created
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others;
  // the "ENZ" suffix stands for "East(X)-North(Y)-Zenith(Z)":
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using DimLessVTopRotENZ  = DimLessV<TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopRotENZ     = PosKV   <TopoCRotCOS_ENZ<BBody, L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopRotENZ     = VelKV   <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopRotENZ      = VelV    <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopRotENZ      = AccV    <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopRotENZ    = ForceV  <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopRotENZ   = AngVelV <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopRot_ENZ  = AngAccV <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopRot_ENZ  = AngMomV <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopRot_ENZ    = TorqV   <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopRot_ENZ     = MoIV    <TopoCRotCOS_ENZ<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopRot_ENZ = MoIRateV<TopoCRotCOS_ENZ<BBody, L>, B>;

  //=========================================================================//
  // "TopoCRotCOS" Struct:                                                   //
  //=========================================================================//
  // Origin   : A point on the Body Surface given by "L";
  // XYZ Axes : Parallel to those of the "BodyCRotCOS" (ie for Earth they are
  //            parallel to the ITRS axes);
  // TimeScale: TT for Earth, TDB otherwise;
  // This COS is HIGHLY NON-INERTIAL: its Origin is rotating (and otherwise mov-
  // ing) in the inertial space, and the axes are rotating as well.
  // Suitable eg for modeling of Launch Ops and for Ground Tracks!
  // XXX: This COS (and its Vectors) do not have any suffix, unlike the above
  // "TopoCRotCOS_ENZ":
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCRotCOS
  {
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    TopoCRotCOS() = delete;
    // No objs of this struct are to be created
  };

  //-----------------------------------------------------------------------//
  // Technical Stuff, for Template Meta-Programming:                       //
  //-----------------------------------------------------------------------//
  // "IsTopoCRotCOS":
  template<typename COS>
  constexpr inline bool IsTopoCRotCOS = false;

  template<Body BBody,  Location<BBody> const* L>
  constexpr inline bool IsTopoCRotCOS<TopoCRotCOS<BBody, L>> = true;

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using DimLessVTopRot = DimLessV<TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopRot    = PosKV   <TopoCRotCOS<BBody,  L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopRot    = VelKV   <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopRot     = VelV    <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopRot     = AccV    <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopRot   = ForceV  <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopRot  = AngVelV <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopRot  = AngAccV <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopRot  = AngMomV <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopRot    = TorqV   <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopRot     = MoIV    <TopoCRotCOS<BBody,  L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopRot = MoIRateV<TopoCRotCOS<BBody,  L>, B>;

  //=========================================================================//
  // Vector Transforms between the "{Body,Topo}CRotCOS"es:                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // BodyC -> TopoC:                                                         //
  //-------------------------------------------------------------------------//
  // Here the CallER needs to specify the "TCOS" explicitly, which must be some
  // "TopoCRotCOS":
  //
  template<typename TCOS, Body  B = Body::UNDEFINED, typename  DQ>
  constexpr Vector3D<DQ,  TCOS, B> ToTopoC
    (Vector3D<DQ, BodyCRotCOS<TCOS::BaseBody>, B> const& a_body_c)
  requires(IsTopoCRotCOS<TCOS>)
  {
    // Position is changed by parallel translation; all other vectors remain un-
    // changed. But this is only possible if both COSes have compatible Time-
    // Stamps:
    TT uniTT = a_body_c.UnifyCOSTSs(TCOS::Loc->PosKV().GetCOSTS());
    return
      IsAnyLen<DQ>
      ? Vector3D<DQ, TCOS, B>
        {
          uniTT,
          a_body_c.x() - TCOS::Loc->PosKV().x(),
          a_body_c.y() - TCOS::Loc->PosKV().y(),
          a_body_c.z() - TCOS::Loc->PosKV().z()
        }
      : Vector3D<DQ, TCOS, B>
        { uniTT, a_body_c.x(), a_body_c.y(), a_body_c.z() };
  }

  //-------------------------------------------------------------------------//
  // TopoC -> BodyC:                                                         //
  //-------------------------------------------------------------------------//
  // Here the "TopoCRotCOS" can be inferred from the argument; the CallER does
  // NOT need to specify any template params:
  //
  template<Body BBody, Location<BBody> const*  L,
           Body B = Body::UNDEFINED, typename  DQ>
     Vector3D<DQ, BodyCRotCOS<BBody>,    B> ToBodyC
    (Vector3D<DQ, TopoCRotCOS<BBody, L>, B> const& a_topo_c)
  {
    // Again, unify the TimeStamps:
    TT uniTT = a_topo_c.UnifyCOSTSs(L->PosKV().GetCOSTS());
    return
      IsAnyLen<DQ>
      ? PosKVRot<BBody, B>
        {
          uniTT,
          L->PosKV().x() + a_topo_c.x(),
          L->PosKV().y() + a_topo_c.y(),
          L->PosKV().z() + a_topo_c.z()
        }
      :
        Vector3D<DQ, BodyCRotCOS<BBody>, B>
        { uniTT, a_topo_c.x(), a_topo_c.y(), a_topo_c.z() };
  }

  //=========================================================================//
  // "TopoCLaunchCOS" Struct:                                                //
  //=========================================================================//
  // Similar to "TopoCRotCOS", but rotated wrt the latter around the Z axis by
  // the angle (counter-clock-wise) (Pi/2-A), where A is the Launch Azimuth.
  // XXX:  Unfortunately, the latter is dynamic (mission-specific),
  // so it cannot be made into a template param:
  // Often, this COS would be the most convenient one for integration of the
  // Equations of Motion during the Launch phase:
  //
  template<Body BBody, Location<BBody> const* L>
  struct TopoCLaunchCOS
  {
    // XXX: This COS is typically Earth-based, but we allow other Bodies as
    // well, so both TT and TDB are supported:
    constexpr static  Body BaseBody              = BBody;
    constexpr static  Location<BBody> const* Loc = L;
    constexpr static  bool HasFixedAxes          = false;
    constexpr static  bool HasFixedOrigin        = false;
    using TimeScale = std::conditional_t<BBody == Body::Earth, TT, TDB>;

    static_assert(L != nullptr);
    // No objs of this struct are to be created:
    TopoCLaunchCOS() = delete;
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors in this COS:                       //
  //-------------------------------------------------------------------------//
  // NB: Using "Len_km" for Pos and Vel Vectors, and "Len_m" for all others:
  //
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using DimLessVTopLaunch = DimLessV<TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using PosKVTopLaunch    = PosKV   <TopoCLaunchCOS<BBody, L>, B>;

  // "VelK" and "Vel":
  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelKVTopLaunch    = VelKV   <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using VelVTopLaunch     = VelV    <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AccVTopLaunch     = AccV    <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using ForceVTopLaunch   = ForceV  <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngVelVTopLaunch  = AngVelV <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngAccVTopLaunch  = AngAccV <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using AngMomVTopLaunch  = AngMomV <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using TorqVTopLaunch    = TorqV   <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIVTopLaunch     = MoIV    <TopoCLaunchCOS<BBody, L>, B>;

  template<Body BBody, Location<BBody> const* L, Body B = Body::UNDEFINED>
  using MoIRateVTopLaunch = MoIRateV<TopoCLaunchCOS<BBody, L>, B>;
}
// End namespave SpaceBallistics
