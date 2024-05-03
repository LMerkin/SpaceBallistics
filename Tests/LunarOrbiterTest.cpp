// vim:ts=2:et
//===========================================================================//
//                             "LunarOrbiter.cpp":                           //
//                  Integration of the Lunar Orbiter Motion                  //
//                    in an Irregualr Gravitational Field                    //
//===========================================================================//
#include "SpaceBallistics/PhysForces/GravityField.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/Locations.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <iostream>

using namespace SpaceBallistics;
using namespace std;

namespace
{
  //=========================================================================//
  // ODE RHS, GSL-Compatible:                                                //
  //=========================================================================//
  // The dimensionality of the ODE system to be solved is 6:
  constexpr static int ODEDim = 6;

  // XXX:
  // (*) For GSL compatibility reasons, the args of this function are NOT
  //     dimensioned; however, all internal computations use DimTypes;
  // (*) Currently, only the (quite complex)  Lunar Gravity Field is used
  //     to compute the RHS; Solar, Earth and Planetary perturbations, as
  //     well as the effects of non-inertiality of the SelenoCentricFixed
  //     CoS, are currently OMITTED:
  //
  int ODERHS
  (
    double       a_t,
    double const a_y    [ODEDim],
    double       a_y_dot[ODEDim],
    void*        // UNUSED
  )
  {
    // Co-Ords and Velocity Components in the "quasi-inertial" SelenoCentric
    // Fixed COS:
    PosVBF<Body::Moon> posF{{Len(a_y[0]), Len(a_y[1]), Len(a_y[2])}}; // m
    Time   t   (a_t);                                                 // sec

    // ("UnTyped") derivatives of those Co-Ords are the corresp Velocities:
    a_y_dot[0] = a_y[3];
    a_y_dot[1] = a_y[4];
    a_y_dot[2] = a_y[5];

    // Now compute the Accelerations. To that end, we need to convert "posF"
    // into the Rotating COS, compute the accelerations there,  and  convert
    // them back into the Fixed COS.
    // We assume that at t=0, the instantaneous ("snap-shot") Rotating COS
    // coincides with the Fixed one; then apply the Moon Rotation Angle  :

    // Siderial Rotation Period:
    constexpr Time PMoon  = To_Time(27.321661_day);

    // Moon Rotation Angle:
    double         MRA    = TwoPi<double> * double(t / PMoon);
    double         cosMRA = Cos(MRA);
    double         sinMRA = Sin(MRA);

    // Co-Ords in the Rotating system via those in the Fixed one:
    PosVBR<Body::Moon> posR
    {{
      cosMRA * posF[0] + sinMRA * posF[1],
      cosMRA * posF[1] - sinMRA * posF[0],
      posF[2]
    }};
    // Acceleration in the Rotating System: Must be cleared first:
    AccVBR<Body::Moon> accR {{Acc(0.0), Acc(0.0), Acc(0.0)}};

    // Now try to compute the actual acceleration:
    // collision with the Lunar surface; 
    try
    {
      GravityField<Body::Moon>::GravAcc(Time(a_t), posR, &accR);
    }
    catch (GravityField<Body::Moon>::ImpactExn const& exn)
    {
      // An exception would mean a likely collison with Lunar surface; stop the
      // integration immediately:
      cout << exn.m_t.Magnitude() << "  " << To_Len_km(exn.m_h) << endl;
      cout << "# LUNAR SURFACE IMPACT NEAR lambda = "
           << To_Angle_deg(exn.m_lambda) << ", phi = "
           << To_Angle_deg(exn.m_phi)    << endl;
      return GSL_EBADFUNC;
    }
    // If OK: Convert "accR"  back into the Fixed COS:
    AccVBF<Body::Moon> accF
    {{
      cosMRA * accR[0] - sinMRA * accR[1],
      sinMRA * accR[0] + cosMRA * accR[1],
      accR[2],
    }};
    // Put them back into the "UnTypes" C array:
    a_y_dot[3] = accF[0].Magnitude();
    a_y_dot[4] = accF[1].Magnitude();
    a_y_dot[5] = accF[2].Magnitude();

    // All Done!
    return 0;
  }
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  // System Definition: Presumably, for an explicit itegration method, no Jacob-
  // ian of the RHS is required. There are no params either:
  gsl_odeiv2_system ODE { ODERHS, nullptr, ODEDim, nullptr };

  // Initial Condition:
  // We assume that at t0=0, the Fixed and Rotating COSes coincide; the Orbiter
  // is in the circular polar orbit around the Moon, over the point (lambda=0,
  // phi=0), at the altitude "h", moving North:
  //
  constexpr Time t0     = 0.0_sec;
  constexpr Len  h0     = To_Len(20.0_km);
  constexpr Len  ReMoon = Location    <Body::Moon>::Re;
  constexpr GM   KMoon  = GravityField<Body::Moon>::K;
  constexpr Len  r0     = ReMoon     + h0;
  constexpr Vel  V0     = SqRt(KMoon / r0);

  // "UnTyped" initial state vector for GSL:
  double y[ODEDim]
  {
    r0.Magnitude(), 0.0, 0.0,
    0.0, 0.0, V0.Magnitude()
  };

  // Run the RKF45 Integrator for 1 Year with 10 sec initial TimeStep:
  constexpr Time   tau     = 10.0_sec;
  constexpr Time   T       = t0 + To_Time(365.25_day);
  // Absolute Precision:
  constexpr Len    AbsPrec = 1.0_m;
  // Relative Precision:
  constexpr double RelPrec = 1e-9;
  // Observation Time Step:
  constexpr Time   tauObs  = 10.0 * tau;

  gsl_odeiv2_driver* ODEDriver =
    gsl_odeiv2_driver_alloc_y_new
      (&ODE,            gsl_odeiv2_step_rkf45,
       tau.Magnitude(), AbsPrec.Magnitude(), RelPrec);
  assert(ODEDriver != nullptr);

  // TIME-MARSHALLING:
  double t = t0.Magnitude();
  while (t < T.Magnitude())
  {
    double t1 =  t + tauObs.Magnitude();
    int    rc =  gsl_odeiv2_driver_apply(ODEDriver, &t, t1, y);

    if (UNLIKELY(rc != 0))
    {
      cout << "# ERROR, exiting..." << endl;
      break;
    }
    // Output the current Altitude:
    Len_km  h =
      To_Len_km(Len(SqRt(Sqr(y[0]) + Sqr(y[1]) + Sqr(y[2]))) - ReMoon);
    cout << t << "  " << h << endl;
  }

  // De-Allocate the Driver:
  (void) gsl_odeiv2_driver_free(ODEDriver);
  return 0;
}
