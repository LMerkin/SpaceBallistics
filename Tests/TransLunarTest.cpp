// vim:ts=2:et
//===========================================================================//
//                         "Tests/TransLunarTest.cpp":                       //
//===========================================================================//
// FIXME: This implementation uses a GPL integrator which is NOT DimTypes-
// aware. Must be replaced by our own DimTypes-based integrator:
//    
#include "SpaceBallistics/CoOrds/KeplerOrbits.hpp"
#include "SpaceBallistics/PhysForces/BodyData.hpp"
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/Maths/Dichotomy.hpp"
#include "SpaceBallistics/Utils.hpp"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <cstdio>

using namespace SpaceBallistics;
using namespace std;

namespace
{
  //=========================================================================//
  // "Context" Struct:                                                       //
  //=========================================================================//
  struct Context
  {
    Time                    m_mjs;    // TimeStamp (MJS_TDB)
    PosKV_GCRS<>            m_pos;    // SC Pos    (DE440T)
    VelKV_GCRS<>            m_vel;    // SC Vel    (DE440T)
    PosKV_GCRS<Body::Moon>  m_posM;   // Moon Pos  (DE440T)
    VelKV_GCRS<Body::Moon>  m_velM;   // Moon Vel  (DE440T)
    LenK                    m_rho;    // SC Distance from Earth
    LenK                    m_rhoM;   // SC Distance from the Moon (DE440T)
  };

  //=========================================================================//
  // "MkTLO":                                                                //
  //=========================================================================//
  // Computes an Elliptical Keplerian Approximation for the TransLunar Traject-
  // ory (Orbit), in GCRS:
  //
  EllipticOrbit<GCRS> MkTLO
  (
    // "a_t0" is the Evaluation Time. The Lunar "rendez-vois" of the SC will
    // occur at the first available opportunity after the last Moon  perigee
    // which occurred not later than "a_t0":
    TDB   a_t0,

    // "a_fQ": Mean Anomaly of the Moon which corresponds to the ApoGee of the
    // TLO; the Lunar "rendez-vois" of the SC will occur near that point. From
    // the Delta-V considerations, it would be advantageous to have "a_fQ" near
    // 0 (near the Moon Orbit perigee), but it may not always be possible  due
    // to various launch constraints:
    //
    Angle a_fQ     = 0.0_rad,

    // "a_deltaQ" is the excess ApoGee distance of the LRO over the Moon Orbit
    // Radius-Vector at "a_fQ"; eg for SelenoCentric Orbit (SCO) Insertion, it
    // should be greater than the Moon Radius:
    //
    LenK  a_deltaQ = 2000.0_km,

    // For SCO insertion, provide some estimated separation from the leading
    // edge of the Moon:
    LenK  a_edge   = 300.0_km,

    // TLO insertion is assumed to be from a circular LEO. TODO: Provide full
    // LEO params here. For the moment, we just provide the curcular altitude
    // relative to the Earth Equator:
    //
    LenK  a_leoH   = 200.0_km
  )
  {
    //-----------------------------------------------------------------------//
    // Instantaneous Keplerian Orbit of the Moon (@ a_t0):                   //
    //-----------------------------------------------------------------------//
    PosKV_GCRS<Body::Moon> posM;
    VelKV_GCRS<Body::Moon> velM;
    DE440T::GetMoonGEqPV(a_t0, &posM, &velM);

    EllipticOrbit<GCRS, Body::Moon> moonOrbit(a_t0, posM, velM);

    //-----------------------------------------------------------------------//
    // Elements of the TLO:                                                  //
    //-----------------------------------------------------------------------//
    // The TLO plane will be the same as the Instantaneous Lunar Orbital Plane:
    Angle  I     = moonOrbit.I();
    Angle  Omega = moonOrbit.Omega();

    // Longitude of the TLO ApoGee (from the common Ascending Node of the Lunar
    // Orbit and the TLO):
    Angle  lQ    = moonOrbit.omega() + a_fQ;

    // Then the Longitude of the TLO perigee from the Node is actually the Arg-
    // ument of the perigee:
    Angle  omega = lQ - PI;

    // Lunar Radius-Vector at "a_fQ" and the TLO ApoGee:
    LenK   Q     = moonOrbit.GetR(a_fQ) + a_deltaQ;

    // XXX: For the moment, we assume that the TLO inserion point is the TLO
    // perigee, corresponding to the LEO radius. Then:
    LenK   q     = BodyData<Body::Earth>::Re + a_leoH;

    // Then the Eccentricity and Semi-Major Axis:
    double e     = double((Q - q)/(Q + q));
    assert(0.0 < e &&   e < 1.0);
    LenK   a     = Q / (e + 1.0);

    // Construct the "provisional" "EllipticOrbit" object for the TLO (with a
    // dummy perigee Time as yet):
    EllipticOrbit<GCRS> tlo0(a, e, I, Omega, omega, TDB{});

    // Finally, the perigee Time. This will be computed relative to the SCO in-
    // sertion time. To determine (approximately) the latter, first compute the
    // intersection point between the Moon Orbit and the TLO  (the one with the
    // SMALLER longitude):
    // "M" is for the Moon, no subscript is for SC:
    //   pM / (1 + eM * cos(l - omegaM)) = p / (1 + e * cos(l - omega)) ;
    // solving this for l=lR, where lR ~< lQ = omega + PI, and therefore:
    //
    LenK   pM     = moonOrbit.p();
    double eM     = moonOrbit.e();
    Angle  omegaM = moonOrbit.omega();
    LenK   p      = tlo0.p();

    Angle  lR     =
      Dichotomy
      (
        [pM, eM, omegaM, p, e, omega](Angle a_l) -> LenK
        {
          return p  * (1.0 + eM * Cos(double(a_l - omegaM))) -
                 pM * (1.0 + e  * Cos(double(a_l - omega )));
        },
        omega,   lQ,  Angle(1e-8)
      );
    // "fRM": The True Anomaly of the Moon at "lR":
    Angle fRM =  lR - omegaM;

    // The corresp distance (must be same for Moon and SC):
    LenK  rRM =  moonOrbit.GetR(fRM);

    // Adjust "fRM" by the angle corresponding to "a_edge" + Moon Radius:
    fRM -= Angle(double((BodyData<Body::Moon>::Re + a_edge) / rRM));

    // Next time (after the previous Moon perigee) the Moon will be at that fRM:
    Time tauM = moonOrbit.TimeSincePeriFocus(fRM).first;
    TDB  tRM  = moonOrbit.T() + tauM;

    // The corresp True Anomaly of the TLO:
    Angle fR  = lR - omega;

    // The TLO flight time from the perigee (XXX: it may change!) to the "fR":
    Time tau  = tlo0.TimeSincePeriFocus(fR).first;
    TDB  T    = tRM - tau;

    // Finally, the full TLO Elements:
    return EllipticOrbit<GCRS>(a, e, I, Omega, omega, T);
  }

  //=========================================================================//
  // RHS for the Exact TLO Integration (GSL-Compatible):                     //
  //=========================================================================//
  // The dimensionality of the ODE system to be solved is 6:
  constexpr static int ODEDim = 6;

  // XXX: For GSL compatibility reasons, the args of this function are NOT
  // dimensioned; however, all internal computations use DimTypes:
  int ODERHS
  (
    double       a_mjs_tdb,
    double const a_y    [ODEDim], // km
    double       a_y_dot[ODEDim], // km/sec
    void*        a_ctx            // Context
  )
  {
    // ("UnTyped") derivatives of the Co-Ords are the corresp Vels:
    a_y_dot[0]   = a_y[3];
    a_y_dot[1]   = a_y[4];
    a_y_dot[2]   = a_y[5];

    // Re-construct the TDB:
    Time mjs_tdb{a_mjs_tdb};
    TDB      tdb{  mjs_tdb};

    // SC (index "K") Pos and Vel in the GCRS (ie relative to Earth):
    PosKV_GCRS<> posEK {LenK(a_y[0]), LenK(a_y[1]), LenK(a_y[2])};
    VelKV_GCRS<> velEK {VelK(a_y[3]), VelK(a_y[4]), VelK(a_y[5])};

    // GeoCentric position and Velocity of the Moon:
    PosKV_GCRS<Body::Moon>    posEM;
    VelKV_GCRS<Body::Moon>    velEM;
    DE440T::GetMoonGEqPV(tdb, &posEM, &velEM);

    // Need BaryCentric (BCRS) positions of the Sun and Earth:
    PosKV_BCRS<Body::Earth>   posBE;
    DE440T::GetPlanetBarEqPV<Body::Earth>(tdb, &posBE, nullptr);

    PosKV_BCRS<Body::Sun>     posBS;
    DE440T::GetPlanetBarEqPV<Body::Sun>  (tdb, &posBS, nullptr);

    // Then the GeoCentric position of the Sun:
    PosKV_GCRS<Body::Sun>     posES  = posBS - posBE;

    // Positions of the SC (index "K") relative to the Moon and Sun (STILL in
    // GCRS, not in a SelenoCentric or HelioCentric system):
    PosKV_GCRS<>              posMK = posEK  - posEM;
    PosKV_GCRS<>              posSK = posEK  - posES;

    // Distances from the SC ("K") to Earth, Moon and Sun:
    LenK rhoEK = LenK(posEK);
    LenK rhoMK = LenK(posMK);
    LenK rhoSK = LenK(posSK);

    LenK rhoEM = LenK(posEM);
    LenK rhoES = LenK(posES);

    // Accelerations of the SC by Earth (the main central force):
    auto accEK = (- DE440T::K<Body::Earth> / Cube(rhoEK)) * posEK;

    // IMPORTANT: However, the GeoCentric perturbation accelerations  of the SC
    // due to the Moon and Sun are NOT just accelerations from those resp bodi-
    // es. Rather,  they are DIFFERENCES between those accelerations of the SC
    // and the accelerations of the Earth caused by the Moon and Sun:
    auto prtMK =
      - DE440T::K<Body::Moon> * (posMK / Cube(rhoMK) + posEM / Cube(rhoEM));

    auto prtSK =
      - DE440T::K<Body::Sun>  * (posSK / Cube(rhoSK) + posES / Cube(rhoES));

    // So the total acceleration is:
    auto acc   = accEK + prtMK + prtSK;

    // Put the Accelerations back into the "UnTyped" C array:
    a_y_dot[3] = acc.x().Magnitude();
    a_y_dot[4] = acc.y().Magnitude();
    a_y_dot[5] = acc.z().Magnitude();

    // Memoise some intermediate data in "ctx":
    Context* ctx = reinterpret_cast<Context*>(a_ctx);
    ctx->m_mjs   = mjs_tdb;
    ctx->m_pos   = posEK;
    ctx->m_vel   = velEK;
    ctx->m_posM  = posEM;
    ctx->m_velM  = velEM;
    ctx->m_rho   = rhoEK;
    ctx->m_rhoM  = rhoMK;

    // All Done!
    return 0;
  }
}
 
//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  TDB t0{TT{UTC{{2025, 1, 1}}}};

  //------------------------------------------------------------------------//
  // Construct the TLO:                                                     //
  //------------------------------------------------------------------------//
  // XXX: Once again, it is currently assumed that the TLO
  // insertion occurs at the perigee; the perigee time is in general NOT "t0"!
  EllipticOrbit<GCRS> tlo = MkTLO(t0);

  printf("TLO:\n\ta=%.3lf\n\tP=%.3lf\n\te=%.9lf\n\tQ=%.3lf\n\tI=%.6lf\n"
         "\tOmega=%.6lf\n\tomega=%.6lf\n\tM0=%.6lf\n\tT=%.3lf\n",
         tlo.a()                    .Magnitude(),
         tlo.P()                    .Magnitude(),
         tlo.e(),
         tlo.Q()                    .Magnitude(),
         To_Angle_deg(tlo.I())      .Magnitude(),
         To_Angle_deg(tlo.Omega())  .Magnitude(),
         To_Angle_deg(tlo.omega())  .Magnitude(),
         To_Angle_deg(tlo.M0())     .Magnitude(),
         tlo.T().GetTimeSinceEpoch().Magnitude());

  // Get the initial Pos and Vel vectors for the TLO (at the perigee, hence
  // f=0):
  PosKV_GCRS<> pos0;
  VelKV_GCRS<> vel0;
  tlo.GetPV(0.0_rad, &pos0, &vel0);

  printf("TLO Insertion: R=%.3lf, V=%.3lf\n",
         LenK(pos0).Magnitude(), VelK(vel0).Magnitude());

  //-------------------------------------------------------------------------//
  // System Definition for GSL-based direct TLO integration:                 //
  //-------------------------------------------------------------------------//
  // Presumably, for an explicit integration method, no Jacobian of the RHS is
  // required:
  Context ctx;
  gsl_odeiv2_system ODE { ODERHS, nullptr, ODEDim, &ctx };

  // "UnTyped" initial state vector for GSL:
  double y[ODEDim]
  {
    pos0.x().Magnitude(), pos0.y().Magnitude(), pos0.z().Magnitude(),
    vel0.x().Magnitude(), vel0.y().Magnitude(), vel0.z().Magnitude()
  };

  // Construct the RKF45 Integrator with 10 sec initial TimeStep (used by the
  // Driver internally):
  constexpr Time   tau     = 1.0_sec;
  // Observation Time Step (between external Driver invocations):
  constexpr Time   tauObs  = 10.0 * tau;

  // Absolute Precision:
  constexpr LenK   AbsPrec = To_Len_km(10.0_m);
  // Relative Precision:
  constexpr double RelPrec = 1e-9;

  gsl_odeiv2_driver* ODEDriver =
    gsl_odeiv2_driver_alloc_y_new
      (&ODE,            gsl_odeiv2_step_rkf45,
       tau.Magnitude(), AbsPrec.Magnitude(), RelPrec);
  assert(ODEDriver != nullptr);

  // Time Range of Integration: 10d are enough (from the TLO insertion at the
  // perigee):
  TDB const tFrom = tlo.T();
  TDB const tTo   = tFrom + To_Time(10.0_day);

  //-------------------------------------------------------------------------//
  // Run the Time-Marshalling:                                               //
  //-------------------------------------------------------------------------//
  double       mjs   = tFrom.GetTimeSinceEpoch().Magnitude();
  double const mjsTo = tTo  .GetTimeSinceEpoch().Magnitude();

  while (mjs < mjsTo)
  {
    double mjs1 = mjs + tauObs.Magnitude();


    int rc = gsl_odeiv2_driver_apply(ODEDriver, &mjs, mjs1, y);

// FIXME:
printf("%.3lf\t%.3lf\t%.3lf\n", ctx.m_mjs.Magnitude(), ctx.m_rho.Magnitude(), ctx.m_rhoM.Magnitude());

    if (UNLIKELY(rc != 0))
    {
      cout << "# ERROR: RC=" << rc << ", exiting..." << endl;
      break;
    }
  }
  // De-Allocate the Driver:
  (void) gsl_odeiv2_driver_free(ODEDriver);

  // All Done!
  return 0;
}
