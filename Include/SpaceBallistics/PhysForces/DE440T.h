// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysForces/DE440T.h":                   //
//                         JPL DE440T Ephemerides                            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include "SpaceBallistics/CoOrds/BaryCentricCOS.h"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include <utility>

namespace SpaceBallistics::DE440T
{
  namespace Bits
  {
    //=======================================================================//
    // Objects Provided by DE440T:                                           //
    //=======================================================================//
    // NB: The nomenclature and enumeration of "Object"s is close to, but NOT
    // exactly the same, as that of "Body"es. Instead of "Body::Earth", there
    // is "Object::EMB" (the Earth-Moon System BaryCenter).   Plus, there are
    // some special "Object"s which are not "Body"s:
    //
    enum class Object: int
    {
      // For the following Objects, DE440T provides "BaryCentricCOS" co-ords:
      Sun             = 0,
      Mercury         = 1,
      Venus           = 2,
      EMB             = 3,
      Mars            = 4,
      Jupiter         = 5,
      Saturn          = 6,
      Uranus          = 7,
      Neptune         = 8,
      Pluto           = 9,
      // For the Moon, DE440T provides "GeoCentricFixedCOS" co-ords:
      Moon            = 10,
      // The rest are not co-ords but some special vals:
      EarthNutations  = 11,    // [d(psi),  d(eps)]
      MoonLibrations  = 12,    // [phi, theta, psi]
      TT_TDB          = 13     // [TT-TDB], sec
    };
    constexpr inline int NObjs = int(Object::TT_TDB) + 1;

    //=======================================================================//
    // Some Mass Consts:                                                     //
    //=======================================================================//
    // Earth/Moon Mass Ratio:
    //
    constexpr inline double EMRat = 81.3005682214972154;

    // Gravitational Consts for all Objects.
    // NB: For Mars, Jupiter, Saturn, Uranus, Neptune and Pluto,  the GMK vals
    // provided include the masses of their moons (that is, similar to the EMB
    // convention):
    //
    template<Object Obj>
    constexpr   GMK K;

    template<>
    constexpr inline GMK K<Object::Sun>     = GMK(132712440041.279419);
    template<>
    constexpr inline GMK K<Object::Mercury> = GMK(       22031.868551);
    template<>
    constexpr inline GMK K<Object::Venus>   = GMK(      324858.592   );
    constexpr inline GMK KEarth             = GMK(      398600.435507);
    template<>
    constexpr inline GMK K<Object::Moon>    = GMK(        4902.800118);
    template<>
    constexpr inline GMK K<Object::EMB>     = K<Object::Moon> + KEarth;
    template<>
    constexpr inline GMK K<Object::Mars>    = GMK(       42828.375816);
    template<>
    constexpr inline GMK K<Object::Jupiter> = GMK(   126712764.1     );
    template<>
    constexpr inline GMK K<Object::Saturn>  = GMK(    37940584.8418  );
    template<>
    constexpr inline GMK K<Object::Uranus>  = GMK(     5794556.4     );
    template<>
    constexpr inline GMK K<Object::Neptune> = GMK(     6836527.10058 );
    template<>
    constexpr inline GMK K<Object::Pluto>   = GMK(         975.5     );

    static_assert(KEarth.ApproxEquals(K<Object::Moon> * EMRat, 1e-10));

    //=======================================================================//
    // Implementation Details:                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Representation Consts:                                                //
    //-----------------------------------------------------------------------//
    // The Temporal Range of Records (in TDB):
    constexpr inline TDB From = TDB(Time_day(2'323'696.5)); // 1649-12-18.0
    constexpr inline TDB To   = TDB(Time_day(2'506'320.5)); // 2149-12-21.0

    // The Number of Records:
    constexpr inline int NR   = 5707;

    // Temporal Span of each Record:
    constexpr Time      RecSpan =  To_Time(32.0_day);
    static_assert(double(NR) * RecSpan == To - From);

    // Size of each record in "double"s:
    constexpr inline int ND   = 1122;

    //-----------------------------------------------------------------------//
    // Chebyshev Coeffs Layout for each Object:                              //
    //-----------------------------------------------------------------------//
    // Dim1 (Outer): The Number of SubPeriods per RecSpan:
    constexpr inline int NSPs[NObjs]
      { 2,  4,   2,  2,  1,  1,  1,  1,  1,  1,  8,  4,  4,  8 };

    // Same, as a templated value:
    template<Object Obj>
    constexpr int NSP = NSPs[int(Obj)];

    // Dim2 (Mid):   The Number of CoOrds:
    constexpr inline int NCOs[NObjs]
      {  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  2,  3,  1 };

    // Same, as a templated value:
    template<Object Obj>
    constexpr inline int NCO = NCOs[int(Obj)];

    // Dim3 (Inner): The Number of Chebyshev Coeffs per CoOrd:
    constexpr int NCCs[NObjs]
      { 11, 14, 10, 13, 11,  8,  7,  6,  6,  6, 13, 10, 10, 13 };

    // Same, as a templated value:
    template<Object Obj>
    constexpr int NCC = NCCs[int(Obj)];

    //-----------------------------------------------------------------------//
    // Types of Chebyshev Coeffs (and CoOrds) for different "Object"s:       //
    //-----------------------------------------------------------------------//
    // In most cases, it is "Len_km", but there are some exceptions. We also
    // provide the types for Time Derivatives ("Dots") of the Coeffs,    and
    // the corresp Array Types:
    //
    template<Object Obj>
    struct CTypes
    {
      using T    = Len_km;
      using D    = decltype(T() / 1.0_sec);
      using ArrT = T[3];        // 3D!
      using ArrD = D[3];        // ditto
    };

    template<>
    struct CTypes<Object::EarthNutations>
    {
      using T    = Angle;
      using D    = AngVel;
      // There are 2 Nutation Angles:
      using ArrT = T[NCO<Object::EarthNutations>];
      using ArrD = D[NCO<Object::EarthNutations>];
    };

    template<>
    struct CTypes<Object::MoonLibrations>
    {
      using T    = Angle;
      using D    = AngVel;
      // There are 3 Libration Angles:
      using ArrT = T[NCO<Object::MoonLibrations>];
      using ArrD = D[NCO<Object::MoonLibrations>];
    };

    template<>
    struct CTypes<Object::TT_TDB>
    {
      using T    = Time;
      using D    = decltype(T() / 1.0_sec); // Actually DimLess!
      // TT_TDB are 1D data:
      using ArrT = T[1];
      using ArrD = D[1];
    };

    // The Short-Cuts:
    template<Object Obj>
    using CT    = CTypes<Obj>::T;

    template<Object Obj>
    using DT    = CTypes<Obj>::D;

    template<Object Obj>
    using ArrCT = CTypes<Obj>::ArrT;

    template<Object Obj>
    using ArrDT = CTypes<Obj>::ArrD;

    //-----------------------------------------------------------------------//
    // Data Record Layout:                                                   //
    //-----------------------------------------------------------------------//
    struct Record
    {
      // Time Span of this Record, JD_TDB:
      Time_day const     m_From;
      Time_day const     m_To;

#     ifdef  MK_DE440T_REC_ENTRY
#     undef  MK_DE440T_REC_ENTRY
#     endif
#     define MK_DE440T_REC_ENTRY(ObjName) \
      /* Chebyshev Coeffs: */   \
      CT<Object::ObjName> const m_##ObjName \
        [NSP<Object::ObjName>]  \
        [NCO<Object::ObjName>]  \
        [NCC<Object::ObjName>];

      // So: Chebyshev Coeffs for all Objects, in the following order (which is
      // NOT the same as enumeration of "Object"s  -- Sun comes after the Moon,
      // not before everyone else):
      MK_DE440T_REC_ENTRY(Mercury)
      MK_DE440T_REC_ENTRY(Venus)
      MK_DE440T_REC_ENTRY(EMB)
      MK_DE440T_REC_ENTRY(Mars)
      MK_DE440T_REC_ENTRY(Jupiter)
      MK_DE440T_REC_ENTRY(Saturn)
      MK_DE440T_REC_ENTRY(Uranus)
      MK_DE440T_REC_ENTRY(Neptune)
      MK_DE440T_REC_ENTRY(Pluto)
      MK_DE440T_REC_ENTRY(Moon)
      MK_DE440T_REC_ENTRY(Sun)
      MK_DE440T_REC_ENTRY(EarthNutations)
      MK_DE440T_REC_ENTRY(MoonLibrations)
      MK_DE440T_REC_ENTRY(TT_TDB)
#     undef MK_DE440T_REC_ENTRY
    };
    // The size of the above "Record" is the same as that of "ND" "double"s:
    static_assert(sizeof(Record) == size_t(ND) * sizeof(double));

    //-----------------------------------------------------------------------//
    // Accessors for the Chebyshev Coeffs:                                   //
    //-----------------------------------------------------------------------//
    // Templated Accessor (implemented below by specialisation):
    // "a_s" is the sub-period index:
    template<Object Obj>
    CT<Obj> const* GetChebyshevCoeffs(Record const* a_rec, int a_s);

#   ifdef  MK_DE440T_GET_COEFFS
#   undef  MK_DE440T_GET_COEFFS
#   endif
#   define MK_DE440T_GET_COEFFS(ObjName)       \
    template<> \
    inline CT<Object::ObjName> const*          \
    GetChebyshevCoeffs<Object::ObjName>(Record const* a_rec, int a_s)     \
    { \
      assert(a_rec != nullptr && 0 <= a_s && a_s < NSP<Object::ObjName>); \
      return &(a_rec->m_##ObjName[a_s][0][0]); \
    }
    MK_DE440T_GET_COEFFS(Sun)
    MK_DE440T_GET_COEFFS(Mercury)
    MK_DE440T_GET_COEFFS(Venus)
    MK_DE440T_GET_COEFFS(EMB)
    MK_DE440T_GET_COEFFS(Mars)
    MK_DE440T_GET_COEFFS(Jupiter)
    MK_DE440T_GET_COEFFS(Saturn)
    MK_DE440T_GET_COEFFS(Uranus)
    MK_DE440T_GET_COEFFS(Neptune)
    MK_DE440T_GET_COEFFS(Pluto)
    MK_DE440T_GET_COEFFS(Moon)
    MK_DE440T_GET_COEFFS(EarthNutations)
    MK_DE440T_GET_COEFFS(MoonLibrations)
    MK_DE440T_GET_COEFFS(TT_TDB)
#   undef  MK_DE440T_GET_COEFFS

    //-----------------------------------------------------------------------//
    // The Actual Data:                                                      //
    //-----------------------------------------------------------------------//
    extern double const Data[NR][ND];

    //-----------------------------------------------------------------------//
    // Internal Utils:                                                       //
    //-----------------------------------------------------------------------//
    // TDB ->  (Record, TimeOffset witin that Record):
    // Returns (NULL, 0.0_sec) if the corresp Record was not found:
    //
    std::pair<Record const*, Time> GetRecord(TDB a_t);

    // (Object, Record, TimeOffset within that Record) -> CoOrds:
    template<Object Obj>
    void GetCoOrds
    (
      Record  const* a_record,
      Time           a_dt,
      ArrCT<Obj>     a_pos,           // Output: CoOrds
      ArrDT<Obj>     a_vel = nullptr  // Output: Dots (if provided)
    );
  }
  // End namespace "Bits"

  //=========================================================================//
  // Mass Constants for Bodies:                                              //
  //=========================================================================//
  // Again, for Mars, Jupiter, Saturn, Uranus, Neptune and Pluto, the GMK vals
  // provided include the masses of the corresp planet and its moons (but this
  // convention does NOT hold for Earth; use Bits::K<Object::EMB>  if you need
  // the GMK value for the Earth-Moon System):
  //
  template<Body BodyName>
  constexpr GMK K;

  template<>
  constexpr inline GMK K<Body::Sun>     = Bits::K<Bits::Object::Sun>;
  template<>
  constexpr inline GMK K<Body::Mercury> = Bits::K<Bits::Object::Mercury>;
  template<>
  constexpr inline GMK K<Body::Venus>   = Bits::K<Bits::Object::Venus>;
  template<>
  constexpr inline GMK K<Body::Earth>   = Bits::KEarth;
  template<>
  constexpr inline GMK K<Body::Mars>    = Bits::K<Bits::Object::Mars>;
  template<>
  constexpr inline GMK K<Body::Jupiter> = Bits::K<Bits::Object::Jupiter>;
  template<>
  constexpr inline GMK K<Body::Saturn>  = Bits::K<Bits::Object::Saturn>;
  template<>
  constexpr inline GMK K<Body::Uranus>  = Bits::K<Bits::Object::Uranus>;
  template<>
  constexpr inline GMK K<Body::Neptune> = Bits::K<Bits::Object::Neptune>;
  template<>
  constexpr inline GMK K<Body::Pluto>   = Bits::K<Bits::Object::Pluto>;
  template<>
  constexpr inline GMK K<Body::Moon>    = Bits::K<Bits::Object::Moon>;

  //=========================================================================//
  // "External" Functions: API to using DE440T:                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // BaryCentric Position and Velocity Vectors for Planets and Sun:          //
  //-------------------------------------------------------------------------//
  // NB: Here the param is "Body", not "Object" (which is internal to "Bits"):
  //
  template<Body BodyName>
  void GetPlanetBPV
  (
    TDB        a_tdb,
    PosKVBary* a_pos,      // Output (Position)
    VelKVBary* a_vel       // Output (Velocity); may be NULL
  );

  void GetPlanetBPV
  (
    Body       a_body,     // Same constraints for "a_obj" as above
    TDB        a_tdb,
    PosKVBary* a_pos,      // Output (Position)
    VelKVBary* a_vel       // Output (Velocity); may be NULL
  );

  // A slightly optimised version for the Sun and all Planets:
  // Since the primary use case of this function is to compute the planetary po-
  // sitions for integration of Minor Solar System Bodies' Orbits,  it produces
  // the co-ords of EMB rather than Earth and Moon separately. The output array
  // is:
  // [0: Sun,    1: Mercury, 2: Venus,   3: EMB, 4: Mars, 5: Jupiter,
  //  6: Saturn, 7: Uranus,  8: Neptune, 9: Pluto]:
  //
  void GetPlanetsBPVs
  (
    TDB        a_tdb,
    PosKVBary  a_poss[10], // Output
    VelKVBary  a_vels[10]  // Output (again, may be NULL)
  );

  //-------------------------------------------------------------------------//
  // GeoCentric Position and Velocity of the Moon:                           //
  //-------------------------------------------------------------------------//
  void GetMoonGPV
  (
    TDB         a_tdb,
    PosKVGeoF*  a_pos,
    VelKVGeoF*  a_vel
  );

  //-------------------------------------------------------------------------//
  // Earth Nutations (Long-Period Only):                                     //
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // Moon Librations:                                                        //
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // TDB <-> TT Conversions:                                                 //
  //-------------------------------------------------------------------------//
  TT  TTofTDB(TDB a_tdb);
  TDB TDBofTT(TT  a_tt);

  //-------------------------------------------------------------------------//
  // "SelfTest" (for Temporal Continuity of Data):                           //
  //-------------------------------------------------------------------------//
  void SelfTest();
}
// End namespace SpaceBallistics::DE440T
