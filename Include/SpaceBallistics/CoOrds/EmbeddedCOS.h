// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/CoOrds/EmbeddedCOS.h":                 //
//       Embedded Co-Ords System (Bound to the Physical Axes of LV/SC)       //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/LVSC/LVSC.h"

namespace SpaceBallistics
{
  //=========================================================================//
  // "EmbeddedCOS" Class:                                                    //
  //=========================================================================//
  // Origin: A fixed physical point               in the LV or SC
  // Axes  : Fixed physical axes (eg of symmetry) of the LV or SC
  // NB    : This type is templatised by the LV/SC Type. It is intended to just
  //         stand for itself (no objs of this class are to be created):
  // TimeScale:
  //         XXX: This is debatable. Typically,  the EmbeddedCOS is used during
  //         the LV Ascent-to-Orbit phase, in which case TT should be used. But
  //         it could also be used for some "deep-space" SC operations, in which
  //         case TDB is perhaps more appropriate.  So we use TT by default but
  //         allow overriding it if really required:
  class TT;
  template<LVSC LVSCKind, typename TS = TT>
  struct EmbeddedCOS
  {
    using TimeScale = TS;

    EmbeddedCOS() = delete;	  // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors and Tensors in this COS:           //
  //-------------------------------------------------------------------------//
  template<LVSC LVSCKind, typename TS = TT>
  using PosVEmb     = PosV    <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using VelVEmb     = VelV    <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using AccVEmb     = AccV    <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using ForceVEmb   = ForceV  <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using AngVelVEmb  = AngVelV <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using AngAccVEmb  = AngAccV <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using AngMomVEmb  = AngMomV <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using TorqVEmb    = TorqV   <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using MoIVEmb     = MoIV    <EmbeddedCOS<LVSCKind, TS>>;

  template<LVSC LVSCKind, typename TS = TT>
  using MoIRateVEmb = MoIRateV<EmbeddedCOS<LVSCKind, TS>>;
}
// End namespace SpaceBallistics
