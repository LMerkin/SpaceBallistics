// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysEffects/DE440T.h":                  //
//                    API to the JPL DE440T Ephemerides                      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/CoOrds/BaryCentricCOSes.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/PhysEffects/BodyData.hpp"   // Has DE440T namespace!
#include "SpaceBallistics/PhysEffects/DE440T-Data.h"

namespace SpaceBallistics::DE440T
{
  //=========================================================================//
  // API to DE440T Ephemerides:                                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // "BaryC{Eq,Ecl}COS" Pos and Vel Vectors for Planets and Sun:             //
  //-------------------------------------------------------------------------//
  // In the BaryC Equatorial COS:
  //
  template<Body B>
  void GetPlanetBarEqPV
  (
    TDB             a_tdb,
    PosKV_BCRS<B>*  a_pos,  // Output (Position)
    VelKV_BCRS<B>*  a_vel   // Output (Velocity); may be NULL
  );

  template<Body B>
  void GetPlanetBarEqPV
  (
    Body            a_body, // Same constraints for "a_obj" as above
    TDB             a_tdb,
    PosKV_BCRS<B>*  a_pos,  // Output (Position)
    VelKV_BCRS<B>*  a_vel   // Output (Velocity); may be NULL
  );

  // In the BaryC Ecliptical COS:
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
    PosKV_BCRS<>    a_poss[10], // Output
    VelKV_BCRS<>    a_vels[10]  // Output (again, may be NULL)
  );

  void GetPlanetsBarEclPVs
  (
    TDB             a_tdb,
    PosKVBarEcl<>   a_poss[10], // Output
    VelKVBarEcl<>   a_vels[10]  // Output (again, may be NULL)
  );

  //-------------------------------------------------------------------------//
  // GeoC Equatorial Position and Velocity of the Moon:                      //
  //-------------------------------------------------------------------------//
  // In the GeoC Equatorial Fixed-Axes (GCRS) COS:
  void GetMoonGEqPV
  (
    TDB                         a_tdb,
    PosKV_GCRS<Body::Moon>*     a_pos,  // Output
    VelKV_GCRS<Body::Moon>*     a_vel   // Output (again, may be NULL)
  );

  // In the GeoC Ecliptical Fixed-Axes COS:
  void GetMoonGEclPV
  (
    TDB                         a_tdb,
    PosKVGeoEclFix<Body::Moon>* a_pos,  // Output
    VelKVGeoEclFix<Body::Moon>* a_vel   // Output (again, may be NULL)
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
    constexpr inline TDB      From(Time_day(2'323'696.5)); // 1649-12-18.0
    constexpr inline TDB      To  (Time_day(2'506'320.5)); // 2149-12-21.0
    constexpr inline Time_jyr FromY = 1650.0_jyr;
    constexpr inline Time_jyr ToY   = 2149.0_jyr;

    // Temporal Span of each Record.  NB: "NR", "ND" are in "DE440T-Data.h":
    constexpr Time      RecSpan =  To_Time(32.0_day);
    static_assert(double(NR) * RecSpan == To - From);

    //-----------------------------------------------------------------------//
    // "TemporalConsistencyTest":                                            //
    //-----------------------------------------------------------------------//
    void TemporalConsistencyTest();
  }
  // End namespace Bits
}
// End namespace SpaceBallistics::DE440T
