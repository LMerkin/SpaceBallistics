// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysForces/DE440T.h":                   //
//                         JPL DE440T Ephemerides                            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/Utils.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"

namespace SpaceBallistics::DE440T
{
  //=========================================================================//
  // Constants:                                                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // Gravitational and Mass Params of Bodies:                                //
  //-------------------------------------------------------------------------//
  template<Body B>
  constexpr GMK K;

  template<> constexpr inline GMK K<Body::Sun>     = GMK(132712440041.279419);
  template<> constexpr inline GMK K<Body::Mercury> = GMK(       22031.868551);
  template<> constexpr inline GMK K<Body::Venus>   = GMK(      324858.592   );
  template<> constexpr inline GMK K<Body::Earth>   = GMK(      398600.435507);
  template<> constexpr inline GMK K<Body::Moon>    = GMK(        4902.800118);
  template<> constexpr inline GMK K<Body::EMB>     = K<Body::Earth>
                                                   + K<Body::Moon>;
  template<> constexpr inline GMK K<Body::Mars>    = GMK(       42828.375816);
  template<> constexpr inline GMK K<Body::Jupiter> = GMK(   126712764.1     );
  template<> constexpr inline GMK K<Body::Saturn>  = GMK(    37940584.8418  );
  template<> constexpr inline GMK K<Body::Uranus>  = GMK(     5794556.4     );
  template<> constexpr inline GMK K<Body::Neptune> = GMK(     6836527.10058 );
  template<> constexpr inline GMK K<Body::PlChB>   = GMK(         975.5     );

  // Earth/Moon Mass Ratio:
  constexpr  inline double EMRat      = 81.3005682214972154;
  static_assert(K<Body::Earth>.ApproxEquals(K<Body::Moon> * EMRat, 1e-10));

  //=========================================================================//
  // API to DE440T Ephemerides:                                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "BaryCentric{Eq,Ecl}COS" Pos and Vel Vectors for Planets and Sun:       //
  //-------------------------------------------------------------------------//
  // In the BaryCentric Equatorial COS:
  //
  template<Body B>
  void GetPlanetBarEqPV
  (
    TDB             a_tdb,
    PosKVBarEq<B>*  a_pos,  // Output (Position)
    VelKVBarEq<B>*  a_vel   // Output (Velocity); may be NULL
  );

  template<Body B>
  void GetPlanetBarEqPV
  (
    Body            a_body, // Same constraints for "a_obj" as above
    TDB             a_tdb,
    PosKVBarEq<B>*  a_pos,  // Output (Position)
    VelKVBarEq<B>*  a_vel   // Output (Velocity); may be NULL
  );

  // In the BaryCentric Ecliptical COS (compatible with JPL Horizons):
  template<Body B>
  void GetPlanetBarEclPV
  (
    TDB             a_tdb,
    PosKVBarEcl<B>* a_pos,  // Output (Position)
    VelKVBarEcl<B>* a_vel   // Output (Velocity); may be NULL
  );

  // XXX: In case of run-time Body selection, the output vectors have the gene-
  // tic (UNDEFINED) Body param:
  void GetPlanetBarEclPV
  (
    Body            a_body, // Same constraints for "a_obj" as above
    TDB             a_tdb,
    PosKVBarEcl<>*  a_pos,  // Output (Position)
    VelKVBarEcl<>*  a_vel   // Output (Velocity); may be NULL
  );

  // A slightly optimised version for the Sun and all Planets:
  // Since the primary use case of this function is to compute the planetary po-
  // sitions for integration of Minor Solar System Bodies' Orbits,  it produces
  // the co-ords of EMB rather than Earth and Moon separately. The output array
  // is:
  // [0: Sun,    1: Mercury, 2: Venus,   3: EMB, 4: Mars, 5: Jupiter,
  //  6: Saturn, 7: Uranus,  8: Neptune, 9: PlChB]:
  //
  // XXX: Since the outputs are arrays of Vectors,   we have to use the generic
  // (UNDEFINED) Body param in all those Vectors. As an alternative, once could
  // provide ptrs to the individually-typed Body's Vectors, but this is cumber-
  // some.
  // Again, both Equatorial (J2000.0) and Ecliptical (J2000.0) versions are av-
  // ailable:
  void GetPlanetsBarEqPVs
  (
    TDB             a_tdb,
    PosKVBarEq<>    a_poss[10], // Output
    VelKVBarEq<>    a_vels[10]  // Output (again, may be NULL)
  );

  void GetPlanetsBarEclPVs
  (
    TDB             a_tdb,
    PosKVBarEcl<>   a_poss[10], // Output
    VelKVBarEcl<>   a_vels[10]  // Output (again, may be NULL)
  );

  //-------------------------------------------------------------------------//
  // GeoCentric Equatorial Position and Velocity of the Moon:                //
  //-------------------------------------------------------------------------//
  // In the GeoCentric Equatorial Fixed-Axes (GCRS) COS:
  void GetMoonGEqPV
  (
    TDB                         a_tdb,
    PosKV_GCRS<Body::Moon>*     a_pos,
    VelKV_GCRS<Body::Moon>*     a_vel
  );

  // In the GeoCentric Ecliptical Fixed-Axes COS (compatible with JPL Horizons):
  void GetMoonGEclPV
  (
    TDB                         a_tdb,
    PosKVGeoEclFix<Body::Moon>* a_pos,
    VelKVGeoEclFix<Body::Moon>* a_vel
  );

  //-------------------------------------------------------------------------//
  // Earth Nutations (Long-Period Only) and Moon Librations:                 //
  //-------------------------------------------------------------------------//
  // The time derivatives of Nutations are considered to be negligible in this
  // case:                                [d(psi),  d(eps)]
  void GetEarthNutations(TDB a_tdb, Angle a_nuts [2]);

  // For Moon Librations, time derivatives may still be needed (but may be NULL
  // as well):                            [phi, theta, psi]  ...
  void GetMoonLibrations(TDB a_tdb, Angle a_librs[3], AngVel a_libr_dots[3]);

  //-------------------------------------------------------------------------//
  // TDB <-> TT Conversions:                                                 //
  //-------------------------------------------------------------------------//
  TT  ToTT (TDB a_tdb);
  TDB ToTDB(TT  a_tt);

  //=========================================================================//
  // Implementation Deatils:                                                 //
  //=========================================================================//
  namespace Bits
  {
    //-----------------------------------------------------------------------//
    // The Temporal Range of Data Provided (in TDB):                         //
    //-----------------------------------------------------------------------//
    constexpr inline TDB From(Time_day(2'323'696.5)); // 1649-12-18.0
    constexpr inline TDB To  (Time_day(2'506'320.5)); // 2149-12-21.0

    // The Number of DE440T Data Records:
    constexpr inline int NR  = 5707;

    // Temporal Span of each Record:
    constexpr Time      RecSpan =  To_Time(32.0_day);
    static_assert(double(NR) * RecSpan == To - From);

    // Size of each record in "double"s:
    constexpr inline int ND  = 1122;

    //-----------------------------------------------------------------------//
    // "TemporalConsistencyTest":                                            //
    //-----------------------------------------------------------------------//
    void TemporalConsistencyTest();

    //-----------------------------------------------------------------------//
    // The Actual Data:                                                      //
    //-----------------------------------------------------------------------//
    extern double const Data[NR][ND];
  }
  // End namespace Bits
}
// End namespace SpaceBallistics::DE440T
