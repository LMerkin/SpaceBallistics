// vim:ts=2:et
//===========================================================================//
//                         "Tests/TransLunarTest.cpp":                       //
//===========================================================================//
#include "SpaceBallistics/CoOrds/TimeScales.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.hpp"
#include "SpaceBallistics/CoOrds/TwoBodyOrbit.hpp"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"
#include "SpaceBallistics/PhysEffects/DE440T.h"
#include "SpaceBallistics/Maths/RKF5.hpp"
#include "SpaceBallistics/Maths/Dichotomy.hpp"
#include <iostream>
#include <cstdio>
#include <utility>

using namespace SpaceBallistics;
using namespace std;

namespace
{
  //=========================================================================//
  // "Context" Struct:                                                       //
  //=========================================================================//
  struct Context
  {
    TDB                     m_tdb;    // TimeStamp (TDB)
    PosKV_GCRS<>            m_posEK;  // SC Pos    (DE440T)
    VelKV_GCRS<>            m_velEK;  // SC Vel    (DE440T)
    PosKV_GCRS<Body::Moon>  m_posEM;  // Moon Pos  (DE440T)
    VelKV_GCRS<Body::Moon>  m_velEM;  // Moon Vel  (DE440T)
    LenK                    m_rhoEK;  // SC Distance from Earth
    LenK                    m_rhoMK;  // SC Distance from the Moon (DE440T)
  };

  //=========================================================================//
  // "MkTLO":                                                                //
  //=========================================================================//
  // Computes an Elliptical 2-Body Approximation for the TransLunar Trajectory
  // (Orbit), in GCRS. Returns the TLO and the Moon Orbit:
  //
  std::pair<TwoBodyOrbit<GCRS>, TwoBodyOrbit<GCRS, Body::Moon>> MkTLO
  (
    // "a_t0" is the Evaluation Time. The Lunar "rendez-vois" of the SC will
    // occur at the first available opportunity after the last Moon  perigee
    // which occurred not later than "a_t0":
    TDB   a_t0,

    // "a_fM": Mean Anomaly of the Moon which corresponds to the apogee of the
    // TLO; the Lunar "rendez-vois" of the SC will occur near that point. From
    // the Delta-V considerations, it would be advantageous to have "a_fM" near
    // 0 (near the Moon Orbit perigee), but it may not always be possible  due
    // to various launch constraints:
    //
    Angle a_fM,

    // "a_deltaQ" is the excess apogee distance of the LRO over the Moon Orbit
    // Radius-Vector at "a_fM"; eg for SelenoCentric Orbit (SCO) Insertion, it
    // should be greater than the Moon Radius:
    //
    LenK  a_deltaQ,

    // For SCO insertion, provide some estimated separation from the leading
    // edge of the Moon:
    LenK  a_edge,

    // TLO insertion is assumed to be from a circular LEO. TODO: Provide full
    // LEO params here. For the moment, we just provide the curcular altitude
    // relative to the Earth Equator:
    //
    LenK  a_leoH
  )
  {
    //-----------------------------------------------------------------------//
    // Instantaneous 2-Body GeoCentric Orbit of the Moon (@ a_t0):           //
    //-----------------------------------------------------------------------//
    PosKV_GCRS<Body::Moon> posM;
    VelKV_GCRS<Body::Moon> velM;
    DE440T::GetMoonGEqPV(a_t0, &posM, &velM);

    TwoBodyOrbit<GCRS, Body::Moon> moonOrbit(a_t0, posM, velM);

    //-----------------------------------------------------------------------//
    // Elements of the TLO:                                                  //
    //-----------------------------------------------------------------------//
    // The TLO plane will be the same as the Instantaneous Lunar Orbital Plane:
    Angle  I     = moonOrbit.I();
    Angle  Omega = moonOrbit.Omega();

    // Longitude of the TLO apogee (from the common Ascending Node of the Lunar
    // Orbit and the TLO):
    Angle  lQ    = moonOrbit.omega() + a_fM;

    // Then the Longitude of the TLO perigee from the Node is actually the Arg-
    // ument of the perigee:
    Angle  omega = lQ - PI;

    // Moon Radius-Vector @ "a_fM", and the TLO apogee:
    LenK   Q     = moonOrbit.GetR(a_fM) + a_deltaQ;

    // XXX: For the moment, we assume that the TLO inserion point is the TLO
    // perigee, corresponding to the LEO radius. Then:
    LenK   q     = BodyData<Body::Earth>::Re + a_leoH;

    // Then the Eccentricity and Semi-Major Axis:
    double e     = double((Q - q)/(Q + q));
    assert(0.0 < e &&   e < 1.0);
    LenK   a     = Q / (e + 1.0);

    // Construct the "provisional" "EllipticOrbit" object for the TLO (with a
    // dummy Perigee Time as yet):
    TwoBodyOrbit<GCRS> tlo0 =
      TwoBodyOrbit<GCRS>::MkEllipticOrbit(a, e, I, Omega, omega, TDB{});

    // Finally, the Perigee Time. This will be computed relative to the SCO in-
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

    // Next time (after the previous Moon perigee) the Moon will be at that fRM:
    Time tauM = moonOrbit.TimeSincePeriFocus(fRM);
    TDB  tRM  = moonOrbit.T() + tauM;

    // The corresp True Anomaly of the TLO:
    Angle fR  = lR - omega;

    // The TLO flight time from the perigee to the "fR":
    Time tau  = tlo0.TimeSincePeriFocus(fR);
    TDB  T    = tRM - tau;

    // Now, start a bit earlier (and arrive at the orbit intersection point
    // earlier) in order to avoid direct collision with the Moon and be able
    // to perform the SCO insertion maneuver.
    // The approx time required by the Moon to cover our "edge" (with the Moon
    // radius added):
    //
    a_edge      += BodyData<Body::Moon>::Re;
    Time advTime = a_edge / VelK(velM);
    T           -= advTime;

    // Finally, the full TLO Elements:
    TwoBodyOrbit<GCRS> tlo1 =
      TwoBodyOrbit<GCRS>::MkEllipticOrbit(a, e, I, Omega, omega, T);

    return std::make_pair(tlo1, moonOrbit);
  }

  //=========================================================================//
  // RHS for the Exact TLO Integration:                                      //
  //=========================================================================//
  // State Vector of the SpaceCraft:
  using StateV = tuple<LenK, LenK, LenK, VelK, VelK, VelK>;

  // The ODE RHS is the derivative of "StateV" wrt Time (TDB):
  using RHSV   = tuple<VelK, VelK, VelK, AccK, AccK, AccK>;

  RHSV ODERHS(StateV const& a_s, TDB a_tdb, Context* a_ctx)
  {
    assert(a_ctx != nullptr);

    // Components of the State Vector:
    LenK x  = get<0>(a_s);
    LenK y  = get<1>(a_s);
    LenK z  = get<2>(a_s);
    VelK Vx = get<3>(a_s);
    VelK Vy = get<4>(a_s);
    VelK Vz = get<5>(a_s);

    // XXX: For GCRS and BCRS, we currently use "UnDef"s TimeStamps;
    // SC (index "K") Pos and Vel in the GCRS (ie relative to Earth):
    PosKV_GCRS<> posEK {TT::UnDef(), TT::UnDef(),  x,  y,  z};

    VelKV_GCRS<> velEK {TT::UnDef(), TT::UnDef(), Vx, Vy, Vz};

    // GeoCentric position and Velocity of the Moon:
    PosKV_GCRS<Body::Moon>       posEM;
    VelKV_GCRS<Body::Moon>       velEM;
    DE440T::GetMoonGEqPV(a_tdb, &posEM, &velEM);

    // Need BaryCentric (BCRS) positions of the Sun and Earth:
    PosKV_BCRS<Body::Earth>   posBE;
    DE440T::GetPlanetBarEqPV<Body::Earth>(a_tdb, &posBE, nullptr);

    PosKV_BCRS<Body::Sun>     posBS;
    DE440T::GetPlanetBarEqPV<Body::Sun>  (a_tdb, &posBS, nullptr);

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

    // Memoise some intermediate data in the Context:
    a_ctx->m_tdb   = a_tdb;
    a_ctx->m_posEK = posEK;
    a_ctx->m_velEK = velEK;
    a_ctx->m_posEM = posEM;
    a_ctx->m_velEM = velEM;
    a_ctx->m_rhoEK = rhoEK;
    a_ctx->m_rhoMK = rhoMK;

    // The result:
    return make_tuple(Vx, Vy, Vz, acc.x(), acc.y(), acc.z());
  }

  //=========================================================================//
  // "PropagateTLO":                                                         //
  //=========================================================================//
  struct TLORes
  {
    LenK  m_deltaQ;      // Excess TLO apogee distance
    LenK  m_edge;        // TLO "aiming edge"
    LenK  m_minDist;     // Minimum Distance to the Moon Surface
    TDB   m_minDistTime; // Time of the Closest Approach
    VelK  m_totalV;      // See the impl...
  };  

  TLORes PropagateTLO
  (
    TDB      a_t0,       // Ref Time (TDB)
    Angle    a_fM,       // "Aim Point": True Anomaly of the Moon Orbit
    LenK     a_deltaQ,   // Excess TLO apogee distance
    LenK     a_edge,     // TLO "aiming edge"
    LenK     a_leoH,     // Circular LEO altitude to start from
    Time     a_tau,
    Time     a_tauMax,
    Time     a_tauObs,   // Observation Time Step
    double   a_relPrec,
    Context* a_ctx,
    bool     a_verbose
  )
  {
    assert(a_ctx != nullptr && IsPos(a_tau) && IsPos(a_tauMax) &&
           IsPos(a_tauObs)  && a_tau <= a_tauMax);

    //-----------------------------------------------------------------------//
    // Construct the TLO:                                                    //
    //-----------------------------------------------------------------------//
    // XXX: Once again, it is currently assumed that the TLO insertion occurs
    // at TLO the perigee. We assume that we are fying into the Moon Orbit pe-
    // rigee as well:
    auto [tlo, mOr] = MkTLO(a_t0, a_fM, a_deltaQ, a_edge, a_leoH);

    if (a_verbose)
    {
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
             tlo.T().GetJD().Magnitude());

      printf("Moon:\n\ta=%.3lf\n\tP=%.3lf\n\te=%.9lf\n\tQ=%.3lf\n\tI=%.6lf\n"
             "\tOmega=%.6lf\n\tomega=%.6lf\n\tM0=%.6lf\n\tT=%.3lf\n",
             mOr.a()                    .Magnitude(),
             mOr.P()                    .Magnitude(),
             mOr.e(),
             mOr.Q()                    .Magnitude(),
             To_Angle_deg(mOr.I())      .Magnitude(),
             To_Angle_deg(mOr.Omega())  .Magnitude(),
             To_Angle_deg(mOr.omega())  .Magnitude(),
             To_Angle_deg(mOr.M0())     .Magnitude(),
             mOr.T().GetJD().Magnitude());
    }
    // Get the initial Pos and Vel vectors for the TLO (at the TLO perigee,
    // hence f=0 here):
    PosKV_GCRS<> pos0;
    VelKV_GCRS<> vel0;
    tlo.GetPV(0.0_rad, &pos0, &vel0);

    if (a_verbose)
      printf("TLO Insertion: R=%.3lf, V=%.3lf\n",
             LenK(pos0).Magnitude(), VelK(vel0).Magnitude());

    //-----------------------------------------------------------------------//
    // Time-Marshaling SetUp:                                                //
    //-----------------------------------------------------------------------//
    // Initial state vector:
    StateV s0 { pos0.x(), pos0.y(), pos0.z(), vel0.x(), vel0.y(), vel0.z() };

    // Time Range of Integration: 10d are enough (from the TLO insertion @ the
    // TLO perigee):
    TDB const tFrom = tlo.T();
    TDB const tTo   = tFrom + To_Time(10.0_day);
    TDB       tObs  = tFrom;  // For verbose output, if enabled

    // The RHS capturing the Context:
    auto rhs =
      [a_ctx](StateV const& a_s, Time a_mjs_tdb) -> RHSV
      { return ODERHS(a_s, TDB{a_mjs_tdb}, a_ctx); };

    // The Call-Back:
    // Min distance to the Moon center:
    LenK         minRhoM    = LenK(+Inf<double>);
    TDB          minRhoMTime;
    VelKV_GCRS<> minRhoMVel;  // Relative velocity at the closest approach

    auto cb  =
      [a_ctx, a_verbose, a_tauObs, &tObs, &minRhoM, &minRhoMTime, &minRhoMVel]
      (StateV*, Time, RHSV const&, StateV const&, Time, RHSV const&) -> bool
      {
        // XXX: The actual lambda args are ignored -- all the data are taken
        // from "a_ctx":
        // SC Velocity relative  to the Moon;
        VelKV_GCRS<> velMK = a_ctx->m_velEK - a_ctx->m_velEM;

        if (a_verbose && a_ctx->m_tdb >= tObs)
        {
          printf("%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
                 a_ctx->m_tdb  .GetJD().Magnitude(),
                 a_ctx->m_rhoEK.Magnitude(),
                 a_ctx->m_rhoMK.Magnitude(),
                 VelK(velMK)   .Magnitude());
          tObs += a_tauObs;
        }

        // Minimum separation detection:
        if (a_ctx->m_rhoMK < minRhoM)
        {
          minRhoM     = a_ctx->m_rhoMK;
          minRhoMTime = a_ctx->m_tdb;
          // The SC velocity relative to the Moon at that time:
          minRhoMVel  = a_ctx->m_velEK - a_ctx->m_velEM;
        }
        // Always continue:
        return true;
      };

    //-----------------------------------------------------------------------//
    // Run the Time-Marshaling:                                              //
    //-----------------------------------------------------------------------//
    // XXX: The Time Range (From..To) params must be of the same type as the
    // "tau" steps, so we cannot use TDB for the former, must use plain "Time":
    Time tFin =
      RKF5(&s0,   tFrom.GetTimeSinceEpoch(), tTo.GetTimeSinceEpoch(), rhs,
           a_tau, a_tauMax,  a_relPrec, &cb, nullptr);

    if (UNLIKELY(tFin != tTo.GetTimeSinceEpoch()))
      throw runtime_error("RKF5 Time-Marshaling FAILED");

    //-----------------------------------------------------------------------//
    // The Result:                                                           //
    //-----------------------------------------------------------------------//
    // We need to estimate the Total DeltaV which is made of:
    // (1) the velocity of the TLO perigee (irrespective to how we get there,
    //     eg from the LEO);
    // (2) the velocity relative to the Moon at the closest approach point;
    // (3) minus the circular velocity at that distance from the Moon:
    //
    VelK circV =
      SqRt(DE440T::K<Body::Moon> / (BodyData<Body::Moon>::Re + minRhoM));
 
    TLORes res
    {
      .m_deltaQ      = a_deltaQ,
      .m_edge        = a_edge,
      .m_minDist     = minRhoM - BodyData<Body::Moon>::Re,
      .m_minDistTime = minRhoMTime,
      .m_totalV      = VelK(vel0) + VelK(minRhoMVel) - circV
    };
    return res;
  }
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  //-------------------------------------------------------------------------//
  // Consts:                                                                 //
  //-------------------------------------------------------------------------//
  TDB const       t0{TT{UTC{{2025, 1, 1}}}};
  constexpr Angle fM         = 0.0_rad;
  constexpr LenK  leoH       = 200.0_km;
  bool            verbose    = false;

  // Ranges for "deltaQ" and "edge":
  constexpr LenK deltaQ_from =  5000.0_km;
  constexpr LenK deltaQ_to   = 30000.0_km;
  constexpr LenK deltaQ_step =   250.0_km;

  constexpr LenK edge_from   =  2000.0_km;
  constexpr LenK edge_to     = 15000.0_km;
  constexpr LenK edge_step   =   250.0_km;

  // Params for the RKF5 Integrator:
  constexpr Time   tau       = 1.0_sec;
  constexpr Time   tauMax    = 10.0 * tau;

  // Observation Time Step (between external Driver invocations):
  constexpr Time   tauObs    = 10.0 * tau;

  // Relative Precision    (roughly corresponds to ~40 m at the Moon distance,
  // which is worse than the DE440T accuracy, but OK for our purposes):
  constexpr double RelPrec   = 1e-7;

  //-------------------------------------------------------------------------//
  // Iterate over the "deltaQ" and "edge":                                   //
  //-------------------------------------------------------------------------//
  Context ctx;

  for (LenK deltaQ = deltaQ_from; deltaQ <= deltaQ_to; deltaQ += deltaQ_step)
  for (LenK edge   = edge_from;   edge   <= edge_to;   edge   += edge_step)
  {
    TLORes res =
      PropagateTLO(t0,   fM, deltaQ, edge, leoH, tau, tauMax, tauObs, RelPrec,
                   &ctx, verbose);

    printf("%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
           deltaQ.Magnitude(), edge.Magnitude(), res.m_minDist.Magnitude(),
           res.m_totalV.Magnitude());
  }
  // All Done!
  return 0;
}
