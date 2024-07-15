// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/CoOrds/BaryCentricCOS.h":                  //
//                    The Inertial Co-Ordinate System                        //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "BaryCentricCOS" Class:                                                 //
  //=========================================================================//
  // Origin: Solar System BaryCenter
  // Axes  : ICRS/ICRF
  //
  class BaryCentricCOS
  {
    BaryCentricCOS() = delete;   // No objects construction at all!
  };

  //-------------------------------------------------------------------------//
  // Position, Velocity and other Vectors and Tensors in this COS:           //
  //-------------------------------------------------------------------------//
  using PosVBary    = PosV  <BaryCentricCOS>;
  using VelVBary    = VelV  <BaryCentricCOS>;
  using AccVBary    = AccV  <BaryCentricCOS>;
  using ForceVBary  = ForceV<BaryCentricCOS>;

  // XXX: Probably no point in considering the MOI Tensors and Rotational Vecs
  // in this COS yet...
}
// End namespace SpaceBallistics
