// vim:ts=2:et
//===========================================================================//
//                       "Tests/LunarOrbiterTest.cpp":                       //
//                  Integration of the Lunar Orbiter Motion                  //
//                    in an Irregualr Gravitational Field                    //
//===========================================================================//
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/CoOrds/Locations.h"
#include "SpaceBallistics/PhysEffects/GravityFld.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"
#include <iostream>

using namespace SpaceBallistics;
using namespace std;

namespace
{
  //=========================================================================//
  // Consts:                                                                 //
  //=========================================================================//
  constexpr LenK ReMoon = Location  <Body::Moon>::Re;
  constexpr GMK  KMoon  = GravityFld<Body::Moon>::K;

  //=========================================================================//
  // ODE RHS:                                                                //
  //=========================================================================//
  // XXX: Currently, only the (quite complex)  Lunar Gravity Field  is used to
  //     compute the RHS; Solar, Earth and Planetary perturbations, as well as
  //     the effects of non-inertiality of the SelenoCFixed COS, are currently
  //     OMITTED.
  // The State Vector is:
  using StateV = tuple<LenK, LenK, LenK, VelK, VelK, VelK>;

  // The ODE RHS is the Time derivative of "StateV":
  using RHSV   = tuple<VelK, VelK, VelK, AccK, AccK, AccK>;

  RHSV ODERHS(StateV const& a_s, Time a_t)
  {
    // Co-Ords and Velocity Components in the "quasi-inertial" SelenoCentric
    // Equatorial Fixed-Axes COS ("EqFixCOS"):
    // FIXME: For the moment, we use UnDef as the TimeStamp of all COSes and
    // Vectors, to avoid conversion of "a_t" back into TDB;  "posF" is a tmp
    // vector anyway:
    //
    // Co-Ords:
    LenK x =  get<0>(a_s);
    LenK y =  get<1>(a_s);
    LenK z =  get<2>(a_s);

    // Velocities:
    VelK Vx = get<3>(a_s);
    VelK Vy = get<4>(a_s);
    VelK Vz = get<5>(a_s);

    PosKVEqFix<Body::Moon> posF{TDB::UnDef(), TDB::UnDef(), x, y, z};

    // Compute the Accelerations. To that end, we need to convert "posF"
    // into the Rotating COS, compute the accelerations there,  and  convert
    // them back into the EqFixCOS.
    // We assume that at t=0, the instantaneous ("snap-shot") Rotating COS
    // coincides with the Fixed one; then apply the Moon Rotation Angle.
    // Siderial Rotation Period:
    constexpr Time PMoon  = To_Time(27.321661_day);

    // Moon Rotation Angle:
    double    MRA    = TwoPi<double> * double(a_t / PMoon);
    double    cosMRA = Cos(MRA);
    double    sinMRA = Sin(MRA);

    // Co-Ords in the Rotating system via those in the Fixed one:
    PosKVRot<Body::Moon> posR
    (
      TDB::UnDef(),
      TDB::UnDef(),
      cosMRA * posF[0] + sinMRA * posF[1],
      cosMRA * posF[1] - sinMRA * posF[0],
      posF[2]
    );
    // Acceleration in the Rotating System: Must be cleared first:
    AccKVRot<Body::Moon> accR
      {TDB::UnDef(), TDB::UnDef(), AccK(0.0), AccK(0.0), AccK(0.0)};

    // Now try to compute the actual acceleration; NB: "ImpactExn" is thrown
    // in case of a collision with the Lunar surface:
    try
    {
      GravityFld<Body::Moon>::GravAcc(a_t, posR, &accR);
    }
    catch (GravityFld<Body::Moon>::ImpactExn const& exn)
    {
      // An exception would mean a likely collison with Lunar surface; stop the
      // integration immediately:
      cout << exn.m_t.Magnitude() << "  " << To_Len_km(exn.m_h) << endl;
      cout << "# LUNAR SURFACE IMPACT NEAR lambda = "
           << To_Angle_deg(exn.m_lambda) << ", phi = "
           << To_Angle_deg(exn.m_phi)    << endl;
      // Re-throw the exception to stop Time Marshaling:
      throw;
    }
    // If OK: Convert "accR"  back into the EqFixCOS:
    AccKVEqFix<Body::Moon> accF
    (
      TDB::UnDef(),
      TDB::UnDef(),
      cosMRA * accR[0] - sinMRA * accR[1],
      sinMRA * accR[0] + cosMRA * accR[1],
      accR[2]
    );

    // The result:
    return make_tuple(Vx, Vy, Vz, accF[0], accF[1], accF[2]);
  }
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  // We assume that at t0=0, the Fixed and Rotating COSes coincide; the Orbiter
  // is in the circular polar orbit around the Moon, over the point (lambda=0,
  // phi=0), at the altitude "h", moving North:
  //
  constexpr LenK h0  = 20.0_km;
  constexpr Time t0  = 0.0_sec;
  constexpr LenK r0  = ReMoon     + h0;
  constexpr VelK V0  = SqRt(KMoon / r0);

  // The Initial State:
  StateV s { r0, 0.0_km, 0.0_km, VelK(0.0), VelK(0.0), V0 };

  // Run the RKF45 Integrator for 1 Year with 10 sec initial TimeStep:
  constexpr Time tau = 1.0_sec;
  constexpr Time T   = t0 + To_Time(365.25_day);

  // Relative Precision (equivalent to ~1 m of absolute precision in R):
  constexpr double RelPrec = 5e-7;

  // ODE Solver Call-Back:
  Time nextOutput = t0;

  auto ODECB =
    [&nextOutput](StateV const& a_s, Time a_t) -> bool
    {
      // Output the current Altitude above the Moon surface -- every 100 sec:
      LenK x = get<0>(a_s);
      LenK y = get<1>(a_s);
      LenK z = get<2>(a_s);

      LenK h = SqRt(Sqr(x) + Sqr(y) + Sqr(z)) - ReMoon;

      if (a_t >= nextOutput)
      {
        cout << a_t.Magnitude() << "  " << h.Magnitude() << endl;
        nextOutput += 100.0_sec;
      }
      // Stop if we have hit the surface:
      return IsPos(h);
    };

  // TIME-MARSHALLING: from t0 to T; "s" is being updated in-place:
  (void) RKF5(&s, t0, T, ODERHS, tau, 10.0*tau, RelPrec, &ODECB, &cerr);

  // All Done!
  return 0;
}
