// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/EmbeddedCOS.h":                 //
//       Embedded Co-Ords System (Bound to the Physical Axes of LV/SC)       //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/LVSC/LVSC.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EmbeddedCOS" Struct:                                                   //
  //=========================================================================//
  // Origin: A fixed physical point               in the LV or SC
  // Axes  : Fixed physical axes (eg of symmetry) of the LV or SC
  // NB    : This type is templatised by the LV/SC Type. It is intended to just
  //         stand for itself (no objs of this struct are to be created):
  // TimeScale:
  //         XXX: This is debatable. Typically,  the EmbeddedCOS is used during
  //         the LV Ascent-to-Orbit phase, in which case TT should be used. But
  //         it could also be used for some "deep-space" SC operations, in which
  //         case TDB is perhaps more appropriate.  So we use TT for now.
  // IMPORTANT: For the avoidance of doubt, "EmbeddedCOS" is a SnapShot of such
  //         a COS, taken at some time instant. Otherwise, the LV/SC velocity,
  //         angular velocity, acceleration, angular acceleration etc,  in the
  //         "EmbeddedCOS" would be indentical 0s. Thus, transformation between
  //         the "EmbeddedCOS" and other COSes would always involve the time in-
  //         stant corresponding to which the SnapShot coprresponds.   Ideally,
  //         that time instant should be made an "EmbeddedCOS" param, but it is
  //         not possible to do so statically in C++,  since the SnapShot  Time-
  //         Stamp is usually known at run-time only. For this reason, the Time-
  //         Stamp is installed in the "Vector3D" instead.
  //
  template<LVSC LVSCKind>
  struct EmbeddedCOS
  {
    constexpr static bool HasFixedAxes   = false;
    constexpr static bool HasFixedOrigin = false;
    using TimeScale                      = TT;

    EmbeddedCOS() = delete;	  // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors and Tensors in this COS:           //
  //-------------------------------------------------------------------------//
  // XXX: Since the Origin is an LVSC here, there is probably no need to provide
  // a template param for the Body which is characterised by those Vectors: it
  // is typically the same LVSC, ie "Body::UNDEFINED":
  //
  template<LVSC LVSCKind>
  using PosVEmb     = PosV    <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using VelVEmb     = VelV    <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using AccVEmb     = AccV    <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using ForceVEmb   = ForceV  <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using AngVelVEmb  = AngVelV <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using AngAccVEmb  = AngAccV <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using AngMomVEmb  = AngMomV <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using TorqVEmb    = TorqV   <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using MoIVEmb     = MoIV    <EmbeddedCOS<LVSCKind>>;

  template<LVSC LVSCKind>
  using MoIRateVEmb = MoIRateV<EmbeddedCOS<LVSCKind>>;
}
// End namespace SpaceBallistics
