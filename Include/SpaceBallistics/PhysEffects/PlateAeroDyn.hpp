// vim:ts=2:et
//===========================================================================//
//              "SpaceBallistics/PhysEffects/PlateAeroDyn.hpp":              //
//            AeroDynamics of a Thin Plate at an Angle of Attack             //
//               in a SubSonic, TranSonic and SuperSonic Flow                //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Vector3D.hpp"
#include "SpaceBallistics/Maths/Dichotomy.hpp"
#include <tuple>

namespace SpaceBallistics
{
  //=========================================================================//
  // "PlateAeroDyn":                                                         //
  //=========================================================================//
  // A stand-alone function which computes the AeroDynamic Force acting on a
  // Thin Plate, in any given COS.
  // Returns: (AeroDynForceVec, MachNumber, FullAeroDynCoeff):
  //
  template<typename COS>
  std::tuple<ForceV<COS>, double, double> PlateAeroDyn
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

    // If the "flow" is 0, there is no aerodynamic force. XXX: The TotalCoeff
    // is then undefined; it is used for testing only, anyway:
    if (flow.IsZero())
      return std::make_tuple(ForceV<COS>(), 0.0, NaN<double>);

    // The Speed of Sound:
    Vel       A     = SqRt(a_gamma * a_p / a_rho);
    assert(IsPos(A));

    // The flow Mach number (notionally, @ +oo distance from the Plate):
    Vel       V     = Vel(flow);
    double    M     = double(V / A);
    double    M2    = Sqr(M);
    assert(IsPos(V) && M > 0.0);

    // The aerodynamic force will be directed towards the normal which makes a
    // SHARP angle with the "flow"; we select the normal "n" accordingly:
    Vel       flowN = flow.DotProd(a_n);
    DimLessV<COS> n = IsPos(flowN) ? a_n : - a_n;
    flowN           = Abs  (flowN);
    assert(flowN <= V);

    // The Angle of Attack (0 .. Pi/2):
    double sinAlpha = double(flowN / V);
    double cosAlpha = SqRt(1.0 - Sqr(sinAlpha));

    Angle  alpha (ASin(sinAlpha));
    assert(!IsNeg(alpha) && alpha <= PI_2);

    // FIXME: We currently do not allow too large Angles of Attack -- our model
    // cannot handle them. So impose the following constraint   (up to 20.0_deg
    // is certainly OK):
    assert(alpha < 0.35_rad);

    // If the Angle of Attack  is 0, the aerodynamic force is 0 in all regimes,
    // and the TotalCoeff is also 0:
    if (IsZero(alpha))
      return std::make_tuple(ForceV<COS>(), M, 0.0);

    //-----------------------------------------------------------------------//
    // Flow Regimes:                                                         //
    //-----------------------------------------------------------------------//
    if (M < 0.8)
    {
      //---------------------------------------------------------------------//
      // SubSonic (but still Compressible) Flow:                             //
      //---------------------------------------------------------------------//
      // NB: The LIFT Force is orthogonal to the "flow", NOT exactly collinear
      // with "n"! There is notionally NO DRAG.
      // The Lift Coeff @ M=0 is 2*Pi*alpha; for M > 0, compressibility effects
      // are treated using the von Karman -- Tsien approximation:
      //
      double cL0 = TwoPi<double> * double(alpha);
      double s   = SqRt(1.0 - M2);
      double cL  = cL0 / (s - 0.5 * M2 * cL0 / (1.0 + s));

      // We should normally have cL > 0;   if not (which may happen for M close
      // enoygh to 1, with sufficiently large "alpha"), then this approximation
      // is not good enough:
      assert(cL > 0.0);

      // The absolute value of the Lift Force:
      Force  L   = cL * a_rho * Sqr(V) / 2.0 * a_S;

      // And the Lift Force Vector, decomposed over the "n" and "V" vectors:
      ForceV<COS> F = L / cosAlpha * (n - (sinAlpha / V) * flow);
      assert(F.EuclidNorm().ApproxEquals(L));
      return std::make_tuple(F, M, cL);
    }
    else
    if (M > 1.2)
    {
      //---------------------------------------------------------------------//
      // SuperSonic Flow:                                                    //
      //---------------------------------------------------------------------//
      // In this case, both LIFT (orthogonal to "flow") and DRAG (collinear with
      // "flow") are present. We actually don't need the separate Lift and Drag
      // components; rather,  we compute the total aerodynamic force vector, in
      // the "n" direction.
      // First, the Oblique Shock Wave above the Plate. Solve the equation for
      // the attached weak shock wave angle, in the range [xL, xR] which is de-
      // fined as follows:
      double a  = (M2 - 1.0) * (M2 * (a_gamma - 1.0) + 2.0);
      double b  = 4.0 +   M2 * (M2 * (a_gamma + 1.0) + 2.0 * (a_gamma - 1.0));
      double c  =               M2 * (a_gamma + 1.0) + 2.0;
      double D  = Sqr(b) + 4.0 * a * c;             // Yes, "+"!

      assert(a > 0.0 && b > 0.0 && c > 0.0 && D > b);
      double xL = 1.0 / SqRt(M2 - 1.0);             // Mach Angle
      double xR = SqRt((b + SqRt(D)) / (2.0 * a));  // ArgMax of f(x) (below)
      assert(xL < xR);

      double y  = 0.5 *     sinAlpha / cosAlpha;    // tan(alpha)/2
      assert(y >  0.0);                             // Because alpha > 0

      // Will solve the equation f(x) = 0, 0 <= x < xM:
      auto f =
        [M2, a_gamma, y](double a_x) -> double
        {
          double x2 = Sqr(a_x);
          return
            ((M2 - 1.0) *  x2 - 1.0) /
            (a_x * (M2  * ((a_gamma - 1.0) * x2 + a_gamma + 1.0) +
                   2.0  * (x2 + 1.0)))
            - y;
        };
      // If "f" does not have different signs at the interval boundaries, there
      // is no solution -- this may happen if "alpha" is too large, and instead
      // of an oblique shock weak wave we have a detached strong shock wave  in
      // front of the Plate. So check for the method applicability:
      //
      assert(f(xL) < 0.0 && f(xR) > 0.0);

      // OK, perform a simple dichotomic search for the solution:
      double x = Dichotomy(f, xL, xR, DefaultTol<double>);

      // x = tan(beta), where "beta" is the obliquity of the shock wave:
      double sinBeta  = x / SqRt(Sqr(x) + 1.0);
      double Mn2      = Sqr(M * sinBeta);

      // The pressure ratio underneath the Plate. NB: We normally have Mn2 > 1
      // and p21 > 1, but this is not so in some margin cases:
      double p21      = 1.0 +  2.0 * a_gamma / (a_gamma + 1.0) * (Mn2 - 1.0);

      // Now the Prandtl-Meyer Expansion Flow above the Plate:
      double e        = SqRt((a_gamma + 1.0) / (a_gamma - 1.0));
      auto   nu       =
        [e](double a_M) -> double
        {
          assert(a_M > 1.0);
          double sM  = SqRt(Sqr(a_M) - 1.0);
          return e   * ATan(sM  / e) - ATan(sM);
        };
      // nu(1) = 0, nu(M > 1) is a monotonically-increasing function.
      // Solve the equation nu(Mup) = nu(M) + alpha, so Mup > Mu:
      // Say M=5 is a reasonable upper bound for the Dichotomic search,
      // since above that, we would have a HyperSonic flow anyway, and
      // this model would not be applicable. But to be on a safe side,
      // set M=10 as the upper bound for the root search:
      //
      double z   = nu(M)  +  double(alpha);
      auto   nu0 = [z, &nu] (double a_M) -> double { return nu(a_M) - z; };
      double Mup = Dichotomy(nu0, M, 10.0,  DefaultTol<double>);
      assert(Mup > M);

      // And then the pressure ratio above the Plate:
      double g1  = (a_gamma - 1.0) / 2.0;
      double gp  =  a_gamma / (a_gamma - 1.0);
      double p31 = std::pow((1.0 + g1 * M2) / (1.0 + g1 * Sqr(Mup)), gp);
      assert(p31 < 1.0);

      // Finally, the total aerodynamic force (in the direction of "n"). It
      // should normally be positive in that direction:
      assert(p31 < p21);
      Force  T   = (p21 - p31) * a_p * a_S;

      // The effective TotalCoeff:
      double cR  = double(T / (a_rho * Sqr(V) / 2.0 * a_S));
      return std::make_tuple(T * n, M, cR);
    }
    else
    {
      //---------------------------------------------------------------------//
      // If we got here: TranSonic Flow:                                     //
      //---------------------------------------------------------------------//
      // XXX: There is no easy-to-apply theory in this case. We will simply
      // provide a polynomial approximation, matching the h
      //
      return std::make_tuple(ForceV<COS>(), M, 0.0);
    }
    __builtin_unreachable();
  }
}
// End namespace SpaceBallistics
