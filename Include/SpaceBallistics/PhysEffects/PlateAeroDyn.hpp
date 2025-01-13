// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/PhysEffects/PlateAeroDyn.hpp":              //
//            AeroDynamics of a Thin Plate at an Angle of Attack             //
//               in a SubSonic, TranSonic and SuperSonic Flow                //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "PlateAeroDyn":                                                         //
  //=========================================================================//
  // A stand-alone function which computes the AeroDynamic Force acting on a
  // Thin Plate, in any given COS:
  //
  template<typename COS>
  ForceV<COS> PlateAeroDyn
  (
    VelV<COS>     const& a_v,     // Plate      Velocity in COS
    VelV<COS>     const& a_w,     // Wind (Air) Velocity in COS
    DimLessV<COS> const& a_n,     // Normal Unit Vector  to the Plate
    Pressure             a_p,     // Atmospheric Pressure
    Density              a_rho,   // Atmosphere  Density
    double               a_gamma, // Atmosphere  cP/cV
    Area                 a_S      // Plate Surface Area
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    assert(IsPos(a_p) && IsPos(a_rho) && IsPos(a_S) && a_gamma > 1.0 &&
           ApproxEqual(double(a_n.EuclidNorm()), 1.0));

    //-----------------------------------------------------------------------//
    // Flow wrt the Plate (in COS):                                          //
    //-----------------------------------------------------------------------//
    VelV<COS> flow  = a_w - a_v;

    // If the "flow" is 0, there is no aerodynamic force:
    if (flow.IsZero())
      return ForceV<COS>();

    // The Speed of Sound:
    Vel       a     = SqRt(a_gamma * a_p / a_rho);
    assert(IsPos(a));

    // The flow Mach number (notionally, @ +oo distance from the Plate):
    Vel       v     = Vel(flow);
    double    M     = double(v / a);
    assert(IsPos(v) && M > 0.0);

    // The aerodynamic force will be directed towards the normal which makes a
    // BLUNT angle with the "flow"; we select the normal "n" accordingly:
    Vel       flowN = flow.DotProd(a_n);
    DimLessV<COS> n = IsNeg(flowN) ? a_n : - a_n;

    // The Angle of Attack (0 .. Pi/2):
    Angle alpha(ASin(std::min(double(Abs(flowN) / v), 1.0)));
    assert(!IsNeg(alpha) && alpha <= PI_2);

    // FIXME: We currently do not allow too large Angles of Attack -- our model
    // cannot handle them. So impose the following constraint:
    assert(alpha < 0.5_rad);

    //-----------------------------------------------------------------------//
    // Flow Regimes:                                                         //
    //-----------------------------------------------------------------------//
    if (M < 0.8)
    {
      //---------------------------------------------------------------------//
      // SubSonic (but still Compressible) Flow:                             //
      //---------------------------------------------------------------------//
    }
    else
    if (M > 1.2)
    {
      //---------------------------------------------------------------------//
      // SuperSonic Flow:                                                    //
      //---------------------------------------------------------------------//
    }
    else
    {
      //---------------------------------------------------------------------//
      // TranSonic Flow:                                                     //
      //---------------------------------------------------------------------//
    }
    return ForceV<COS>();
  }
}
// End namespace SpaceBallistics
