// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/PhysEffects/DE440T.hpp":                 //
//                    Functions (but not Data) for DE440T                    //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysEffects/DE440T.h"
#include "SpaceBallistics/Maths/Chebyshev.hpp"
#include "SpaceBallistics/CoOrds/EqEcl.hpp"
#include <cmath>
#include <utility>
#include <cassert>

namespace SpaceBallistics::DE440T
{
  namespace Bits
  {
    //=======================================================================//
    // Objects Provided by DE440T:                                           //
    //=======================================================================//
    // NB: The nomenclature and enumeration of "Object"s is the same as that of
    // "Body"s, with addition of Nutations, Librations and TT_TDB difference,
    // and w/o Earth which is not represented in DE440T directly:
    //
    enum class Object: int
    {
      Sun             = int(Body::Sun),
      Mercury         = int(Body::Mercury),
      Venus           = int(Body::Venus),
      EMB             = int(Body::Earth),    // BEWARE: EMB in place of Earth!
      Mars            = int(Body::Mars),
      Jupiter         = int(Body::Jupiter),
      Saturn          = int(Body::Saturn),
      Uranus          = int(Body::Uranus),
      Neptune         = int(Body::Neptune),
      PlChB           = int(Body::PlChB),    // Pluto-Charon BaryCenter
      Moon            = int(Body::Moon),
      EarthNutations  = int(Body::Moon) + 1, // [d(psi),  d(eps)]
      MoonLibrations  = int(Body::Moon) + 2, // [phi, theta, psi]
      TT_TDB          = int(Body::Moon) + 3  // [TT-TDB], sec
    };
    constexpr inline int NObjs = int(Object::TT_TDB) + 1;

    //-----------------------------------------------------------------------//
    // "Object" from a "Body":                                               //
    //-----------------------------------------------------------------------//
    constexpr Object ObjectOfBody(Body B)
    {
      // NB: Body::Earth is not directly translatable!
      assert(Body::Sun <= B && B <= Body::EMB && B != Body::Earth);
      return
        (B != Body::EMB)
        ? // Then a direct re-interpretation is possible, see "Object" above:
          Object(int(B))
        : Object::EMB;
    }

    //=======================================================================//
    // Implementation Details:                                               //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Chebyshev Coeffs Layout for each Object:                              //
    //-----------------------------------------------------------------------//
    // XXX: IMPORTANT: Earth is not represented directly in DE440T data, so we
    // put 0s in the corresp positions below:
    //
    // Dim1 (Outer): The Number of SubPeriods per RecSpan:
    template<Object Obj>
    int NSP;

    template<> constexpr inline int NSP<Object::Sun>            = 2;
    template<> constexpr inline int NSP<Object::Mercury>        = 4;
    template<> constexpr inline int NSP<Object::Venus>          = 2;
    template<> constexpr inline int NSP<Object::EMB>            = 2;
    template<> constexpr inline int NSP<Object::Mars>           = 1;
    template<> constexpr inline int NSP<Object::Jupiter>        = 1;
    template<> constexpr inline int NSP<Object::Saturn>         = 1;
    template<> constexpr inline int NSP<Object::Uranus>         = 1;
    template<> constexpr inline int NSP<Object::Neptune>        = 1;
    template<> constexpr inline int NSP<Object::PlChB>          = 1;
    template<> constexpr inline int NSP<Object::Moon>           = 8;
    template<> constexpr inline int NSP<Object::EarthNutations> = 4;
    template<> constexpr inline int NSP<Object::MoonLibrations> = 4;
    template<> constexpr inline int NSP<Object::TT_TDB>         = 8;

    // Dim2 (Mid):   The Number of CoOrds. It is 3 in most cases, with the
    // exceptions listed below:
    template<Object Obj>
    constexpr inline int NCO = 3;

    template<> constexpr inline int NCO<Object::EarthNutations> = 2;
    template<> constexpr inline int NCO<Object::TT_TDB>         = 1;

    // Dim3 (Inner): The Number of Chebyshev Coeffs per CoOrd:
    template<Object Obj>
    int NCC;

    template<> constexpr inline int NCC<Object::Sun>            = 11;
    template<> constexpr inline int NCC<Object::Mercury>        = 14;
    template<> constexpr inline int NCC<Object::Venus>          = 10;
    template<> constexpr inline int NCC<Object::EMB>            = 13;
    template<> constexpr inline int NCC<Object::Mars>           = 11;
    template<> constexpr inline int NCC<Object::Jupiter>        =  8;
    template<> constexpr inline int NCC<Object::Saturn>         =  7;
    template<> constexpr inline int NCC<Object::Uranus>         =  6;
    template<> constexpr inline int NCC<Object::Neptune>        =  6;
    template<> constexpr inline int NCC<Object::PlChB>          =  6;
    template<> constexpr inline int NCC<Object::Moon>           = 13;
    template<> constexpr inline int NCC<Object::EarthNutations> = 10;
    template<> constexpr inline int NCC<Object::MoonLibrations> = 10;
    template<> constexpr inline int NCC<Object::TT_TDB>         = 13;

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
      // There are 2 Nutation Angles:  [d(psi), d(eps)]:
      static_assert(NCO<Object::EarthNutations> == 2);
      using ArrT = T[2];
      using ArrD = D[2];
    };

    template<>
    struct CTypes<Object::MoonLibrations>
    {
      using T    = Angle;
      using D    = AngVel;
      // There are 3 Libration Angles: [phi, theta, psi]:
      static_assert(NCO<Object::MoonLibrations> == 3);
      using ArrT = T[3];
      using ArrD = D[3];
    };

    template<>
    struct CTypes<Object::TT_TDB>
    {
      using T    = Time;
      using D    = decltype(T() / 1.0_sec); // Actually DimLess!
      // TT_TDB are 1D data:
      static_assert(NCO<Object::TT_TDB> == 1);
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
      // not before everyone else). HERE THE ORDER IS IMPORTANT!
      MK_DE440T_REC_ENTRY(Mercury)
      MK_DE440T_REC_ENTRY(Venus)
      MK_DE440T_REC_ENTRY(EMB)
      MK_DE440T_REC_ENTRY(Mars)
      MK_DE440T_REC_ENTRY(Jupiter)
      MK_DE440T_REC_ENTRY(Saturn)
      MK_DE440T_REC_ENTRY(Uranus)
      MK_DE440T_REC_ENTRY(Neptune)
      MK_DE440T_REC_ENTRY(PlChB)
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
    MK_DE440T_GET_COEFFS(PlChB)
    MK_DE440T_GET_COEFFS(Moon)
    MK_DE440T_GET_COEFFS(EarthNutations)
    MK_DE440T_GET_COEFFS(MoonLibrations)
    MK_DE440T_GET_COEFFS(TT_TDB)
#   undef  MK_DE440T_GET_COEFFS

    //=======================================================================//
    // "GetRecord":                                                          //
    //=======================================================================//
    // Internal Util. Returns the pair:
    //   (Record Ptr corresponding to the given TDB instant,
    //    Time Offset within the Record):
    //
    inline std::pair<Record const*, Time> GetRecord(TDB a_tdb)
    {
      // XXX: Protection against out-of-range "a_tdb" (this may sometimes happen
      // due to rounding errors):
      assert(From <= a_tdb && a_tdb <= To);

      Time   off = a_tdb - From;
      assert(!IsNeg(off));

      double r0  = double(off / RecSpan);
      double rf  = std::floor(r0);
      int    r   = int(rf);
      Time   dt  = off - rf * RecSpan;

      // The following must hold:
      assert(0 <= r && r < NR && !IsNeg(dt) && dt < RecSpan);

      // We can now extract the Record:
      Record const* rec = reinterpret_cast<Record const*>(&(Data[r][0]));

      // XXX: "rec" should indeed contain the required TDB period, BUT it might
      // happen (due to rounding error in the above computations) that it is off
      // by 1 record:
      //
      TDB from(rec->m_From);
      TDB to  (rec->m_To);
      if (UNLIKELY(a_tdb < from))
      {
        // Move to the left:
        assert(r > 0);
        --rec;
        from = TDB(rec->m_From);
        to   = TDB(rec->m_To);
      }
      else
      if (UNLIKELY(a_tdb > to))
      {
        // Move to the right:
        assert(r < NR-1);
        ++rec;
        from = TDB(rec->m_From);
        to   = TDB(rec->m_To);
      }
      // Re-check the "a_tdb" range:
      assert(from <= a_tdb && a_tdb < to);

      // In any case, re-calculate "dt" for higher precision:
      dt = a_tdb - from;
      assert(!IsNeg(dt) && dt < RecSpan);

      // The result:
      return std::make_pair(rec, dt);
    }

    //=======================================================================//
    // "GetCoOrds":                                                          //
    //=======================================================================//
    // XXX: For internal use only.   The output "a_pos" and "a_vel" are arrays,
    // not vectors, so this function is also suitable for Nutaions, Librations
    // and TT_TDB diffs:
    //
    template<Object Obj>
    void GetCoOrds
    (
      Record const* a_record,           // From a prev call to "GetRecord"
      Time          a_dt,               // ditto
      ArrCT<Obj>    a_pos,              // Non-NULL
      ArrDT<Obj>    a_vel               // May be NULL
    )
    {
      assert(a_record != nullptr && !IsNeg(a_dt) && a_dt < RecSpan &&
             a_pos    != nullptr);

      // Initially "a_dt" is the temporal offset within the whole RecSpan; find
      // the corresp SubPeriod Index ("s1"):
      constexpr Time SP  = RecSpan     / double(NSP<Obj>);  // SubPeriod Length
      double         s0  = double(a_dt / SP);
      double         sf  = std::floor(s0);
      int            s   = int       (sf);
      assert( 0 <=   s   &&  s    <  NSP<Obj>);

      // Compute the offset "tau" from the beginning of the SubPeriod:
      Time   tau  = a_dt - sf * SP;
      assert(tau <= a_dt);

      // In any case, we must get valid "s" and "tau":
      assert(0 <= s && s < NSP<Obj> && !IsNeg(tau) && tau < SP && tau <= a_dt);

      // Compute the argument "x" of Chebyshev Polynomials, in [-1..1]:
      double  x = 2.0 * double(tau / SP) - 1.0;
      assert(-1.0 <= x && x <= 1.0);

      // In case we need the Dots, compute dx/dt:
      auto dxdt = 2.0 / SP;

      // Get the Chebyshev Coeffs for this Object and SubPeriod:
      CT<Obj> const* coeffsBase = GetChebyshevCoeffs<Obj>(a_record, s);

      // For all CoOrds and same "x", evaluate the Chebyshev Sum for the orders
      // [0 .. NCC<Obj>-1]:
      for (int i = 0; i < NCO<Obj>; ++i)
      {
        CT<Obj> const* coeffsI = coeffsBase + i * NCC<Obj>;

        // XXX: BEWARE: The function "Chebyshev::Sum1T" takes the 0th coeff with
        //      factor 0.5, whereas DE440T expansions assume factor 1.0, so this
        //      is corrected explicitly. For "Chebyshev::SumDT", there is no di-
        //      screpancy:
        a_pos[i] = Chebyshev::Sum1T<double, CT<Obj>>(NCC<Obj>-1, coeffsI, x)
                 + 0.5 * coeffsI[0];

        if (a_vel != nullptr)
          // Provide the Dots (Velocities) as well. They are NOT stored  in the
          // DE440T data sets, so we have to differentiate the Chebyshev Series:
          a_vel[i] = Chebyshev::SumDT<double, CT<Obj>>(NCC<Obj>-1, coeffsI, x)
                   * dxdt;
      }
    }

    //=======================================================================//
    // "TemporalConsistencyTest":                                            //
    //=======================================================================//
    void TemporalConsistencyTest()
    {
      TDB expFrom = From;
      TDB expTo;    // Empty as yet

      for (int r = 0; r < NR; ++r)
      {
        Record const* rec = reinterpret_cast<Record const*>(&(Data[r][0]));

        TDB from(rec->m_From);
        TDB to  (rec->m_To);
        expTo  = from + RecSpan;

        assert(from == expFrom && to == expTo);
        expFrom = to;
      }
      assert(expTo == To);
    }
  }
  // End namespace Bits

  //=========================================================================//
  // "GetPlanet[s]BarEqPV[s]":                                               //
  //=========================================================================//
  // In the "BaryCEqCOS", ie in the ICRS/BCRS axes:
  //-------------------------------------------------------------------------//
  // With Compile-Time Object Selection:                                     //
  //-------------------------------------------------------------------------//
  template<Body B>
  void GetPlanetBarEqPV
  (
    TDB             a_tdb,
    PosKV_BCRS<B>*  a_pos,  // Output (Position)
    VelKV_BCRS<B>*  a_vel   // Output (Velocity); may be NULL
  )
  {
    // This is the generic case for the Sun and Planets except Earth (which re-
    // uires a specialised implementations); Moon is also not supported here as
    // we seldom need its BaryC PV; however, EMB is OK here:
    //
    static_assert(B != Body::Earth && B != Body::Moon);
    assert(a_pos != nullptr);

    // Get the record data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    constexpr     Bits::Object  Obj  = Bits::ObjectOfBody(B);
    static_assert(Bits::Object::Sun <= Obj && Obj <= Bits::Object::PlChB);

    // Invoke the generic "GetCoOrds". Due to the DE440T convention, it will
    // yield the PV in the BaryC Eq co-ords directly:
    Bits::GetCoOrds<Obj>
    (
      rec,
      dt,
      a_pos->GetArr(),
      a_vel != nullptr ? a_vel->GetArr() : nullptr
    );
  }

  //-------------------------------------------------------------------------//
  // Specialisation for Earth:                                               //
  //-------------------------------------------------------------------------//
  template<>
  void GetPlanetBarEqPV<Body::Earth>
  (
    TDB                      a_tdb,
    PosKV_BCRS<Body::Earth>* a_pos,  // Output (Position)
    VelKV_BCRS<Body::Earth>* a_vel   // Output (Velocity); may be NULL
  )
  {
    assert(a_pos != nullptr);

    // Get the record data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    // First, get the PV of the Earth-Moon System BaryCenter, in the BaryC Eq
    // System (according to the DE440T convention):
    LenK posEMB[3];
    VelK velEMB[3];
    Bits::GetCoOrds<Bits::Object::EMB>
    (
      rec,
      dt,
      posEMB,
      (a_vel != nullptr) ? velEMB  : nullptr
    );

    // Now get the GeoC PV of the Moon, in the GeoC Eq system (again, according
    // to the DE440T convention):
    LenK posM[3];
    VelK velM[3];
    Bits::GetCoOrds<Bits::Object::Moon>
    (
      rec,
      dt,
      posM,
      (a_vel != nullptr) ? velM : nullptr
    );

    // The GeoC PV of the EMB: Proportional to those of the Moon:
    constexpr double mu = 1.0 / (EMRat + 1.0);

    // So the BaryC Eq PV of Earth will be:
    a_pos->x() = posEMB[0] - mu * posM[0];
    a_pos->y() = posEMB[1] - mu * posM[1];
    a_pos->z() = posEMB[2] - mu * posM[2];

    if (a_vel != nullptr)
    {
      a_vel->x() = velEMB[0] - mu * velM[0];
      a_vel->y() = velEMB[1] - mu * velM[1];
      a_vel->z() = velEMB[2] - mu * velM[2];
    }
  }

  //-------------------------------------------------------------------------//
  // With Run-Time Body Selection:                                           //
  //-------------------------------------------------------------------------//
  // This is just a mutiplexor of the above version, The output verctors have
  // the generic Body::UNDEFINED param:
  //
  void GetPlanetBarEqPV
  (
    Body          a_body,
    TDB           a_tdb,
    PosKV_BCRS<>* a_pos,    // Output (Position)
    VelKV_BCRS<>* a_vel     // Output (Velocity); may be NULL
  )
  {
    switch (a_body)
    {
#   ifdef  DE440T_BODY_CASE
#   undef  DE440T_BODY_CASE
#   endif
#   define DE440T_BODY_CASE(B) \
      case Body::B: \
           /* NB: Cast UnDef Vector ptrs into Body-Specific ones: */ \
           GetPlanetBarEqPV<Body::B> \
             (a_tdb, ToSpecBody<Body::B>(a_pos), ToSpecBody<Body::B>(a_vel)); \
           break;
      DE440T_BODY_CASE(Sun)
      DE440T_BODY_CASE(Mercury)
      DE440T_BODY_CASE(Venus)
      DE440T_BODY_CASE(Earth)
      DE440T_BODY_CASE(Mars)
      DE440T_BODY_CASE(Jupiter)
      DE440T_BODY_CASE(Saturn)
      DE440T_BODY_CASE(Uranus)
      DE440T_BODY_CASE(Neptune)
      DE440T_BODY_CASE(PlChB)
      DE440T_BODY_CASE(EMB)
#   undef  DE440T_BODY_CASE
    default:
      assert(false);
    }
  }

  //-------------------------------------------------------------------------//
  // A slightly optimised version for the Sun and all Major Planets:         //
  //-------------------------------------------------------------------------//
  // Since the primary use case of this function is to compute the planetary po-
  // sitions for integration of Minor Solar System Bodies' Orbits,  it produces
  // the co-ords of EMB rather than Earth and Moon separately. The output array
  // (of Generic vectors) is:
  // [0: Sun,    1: Mercury, 2: Venus,   3: EMB,  4: Mars, 5: Jupiter,
  //  6: Saturn, 7: Uranus,  8: Neptune, 9: PlChB]:
  //
  void GetPlanetsBarEqPVs
  (
    TDB           a_tdb,
    PosKV_BCRS<>  a_poss[10], // Output
    VelKV_BCRS<>  a_vels[10]  // Output (again, may be NULL)
  )
  {
    // Optimisation: The Record is fetched JUST ONCE for all Objects:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

#   ifdef  DE440T_OBJ_PV
#   undef  DE440T_OBJ_PV
#   endif
#   define DE440T_OBJ_PV(I, ObjName) \
    { \
      static_assert(I == int(Bits::Object::ObjName)); \
      Bits::GetCoOrds<Bits::Object::ObjName> \
        (rec, \
         dt,  \
         a_poss  [I].GetArr(), \
         ( a_vels != nullptr ) \
         ? a_vels[I].GetArr()  \
         : nullptr \
      ); \
    }
    DE440T_OBJ_PV(0, Sun)
    DE440T_OBJ_PV(1, Mercury)
    DE440T_OBJ_PV(2, Venus)
    DE440T_OBJ_PV(3, EMB)
    DE440T_OBJ_PV(4, Mars)
    DE440T_OBJ_PV(5, Jupiter)
    DE440T_OBJ_PV(6, Saturn)
    DE440T_OBJ_PV(7, Uranus)
    DE440T_OBJ_PV(8, Neptune)
    DE440T_OBJ_PV(9, PlChB)
#   undef DE440T_OBJ_PV
  }

  //=========================================================================//
  // "GetPlanet[s]BarEclPV[s]":                                              //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // With Compile-Time Body Selection:                                       //
  //-------------------------------------------------------------------------//
  template<Body B>
  void GetPlanetBarEclPV
  (
    TDB             a_tdb,
    PosKVBarEcl<B>* a_pos,    // Output (Position)
    VelKVBarEcl<B>* a_vel     // Output (Velocity); may be NULL
  )
  {
    // Not for the Moon;  but Sun, Planets and EMB are OK:
    static_assert(Body::Sun <= B && B <= Body::EMB && B != Body::Moon);

    // First, compute the PV in the BaryC Eq Co-Ords:
    PosKV_BCRS<B>    posEq;
    VelKV_BCRS<B>    velEq;
    GetPlanetBarEqPV<B>
      (a_tdb, &posEq, a_vel != nullptr ? &velEq : nullptr);

    // Then convert the results into the BaryC Ecl Co-Ords:
    *a_pos = ToEcl<BaryCEclCOS>(posEq);
    if (a_vel != nullptr)
      *a_vel = ToEcl<BaryCEclCOS>(velEq);
  }

  //-------------------------------------------------------------------------//
  // With Run-Time Body selection:                                           //
  //-------------------------------------------------------------------------//
  void GetPlanetBarEclPV
  (
    Body            a_body,   // Same constraints for "a_obj" as above
    TDB             a_tdb,
    PosKVBarEcl<>*  a_pos,    // Output (Position)
    VelKVBarEcl<>*  a_vel     // Output (Velocity); may be NULL
  )
  {
    // First, compute the PV in the BaryC Eq Co-Ords:
    PosKV_BCRS<>    posEq;
    VelKV_BCRS<>    velEq;
    GetPlanetBarEqPV
      (a_body, a_tdb, &posEq, a_vel != nullptr ? &velEq : nullptr);

    // Then convert the results into the BaryC Ecl Co-Ords:
    *a_pos = ToEcl<BaryCEclCOS>(posEq);
    if (a_vel != nullptr)
      *a_vel = ToEcl<BaryCEclCOS>(velEq);
  }

  //-------------------------------------------------------------------------//
  // A slightly optimised version for the Sun and all Major Planets:         //
  //-------------------------------------------------------------------------//
  // Similar to "GetPlanetsBarEqPVs", but for the BaryC Ecl PVs. The output
  // arrays (of Generic vectors) are for:
  // [0: Sun,    1: Mercury, 2: Venus,   3: EMB, 4: Mars, 5: Jupiter,
  //  6: Saturn, 7: Uranus,  8: Neptune, 9: PlChB]:
  //
  void GetPlanetsBarEclPVs
  (
    TDB             a_tdb,
    PosKVBarEcl<>   a_poss[10], // Output
    VelKVBarEcl<>   a_vels[10]  // Output (again, may be NULL)
  )
  {
    // First, get the PVs in the BaryC Eq Co-Ords:
    PosKV_BCRS<>    possEq[10];
    VelKV_BCRS<>    velsEq[10];
    GetPlanetsBarEqPVs(a_tdb, possEq, (a_vels != nullptr) ? velsEq : nullptr);

    // Then convert the results into the BaryC Ecl Co-Ords:
    for (int i = 0; i < 10; ++i)
    {
      a_poss[i] = ToEcl<BaryCEclCOS>(possEq[i]);
      if (a_vels != nullptr)
        a_vels[i] = ToEcl<BaryCEclCOS>(velsEq[i]);
    }
  }

  //=========================================================================//
  // "GetMoonG{Eq,Ecl}PV":                                                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // In the GeoC Eq Fixed-Axes (ICRS/GCRS) COS:                              //
  //-------------------------------------------------------------------------//
  void GetMoonGEqPV
  (
    TDB                     a_tdb,
    PosKV_GCRS<Body::Moon>* a_pos,
    VelKV_GCRS<Body::Moon>* a_vel
  )
  {
    assert(a_pos != nullptr);

    // Get the record data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    // Invoke the generic "GetCoOrds". Due to the DE440T convention, it will
    // yield the PV in the GeoC Eq Fixed-Axes (ICRS/GCRS) coords directly:
    Bits::GetCoOrds<Bits::Object::Moon>
    (
      rec,
      dt,
      a_pos->GetArr(),
      a_vel != nullptr ? a_vel->GetArr() : nullptr
    );
  }

  //-------------------------------------------------------------------------//
  // In the GeoC Ecl Fixed-Axes COS:                                         //
  //-------------------------------------------------------------------------//
  void GetMoonGEclPV
  (
    TDB                         a_tdb,
    PosKVGeoEclFix<Body::Moon>* a_pos,
    VelKVGeoEclFix<Body::Moon>* a_vel
  )
  {
    // First, compute the PV in the GeoC Eq Fixed-Axes (ICRS/GCRS) Co-Ords:
    PosKV_GCRS<Body::Moon> posEq;
    VelKV_GCRS<Body::Moon> velEq;

    GetMoonGEqPV(a_tdb, &posEq, a_vel != nullptr ? &velEq : nullptr);

    // Then convert the results into the GeoC Ecl Co-Ords:
    *a_pos = ToEcl<GeoCEclFixCOS>(posEq);
    if (a_vel != nullptr)
      *a_vel = ToEcl<GeoCEclFixCOS>(velEq);
  }

  //=========================================================================//
  // "ToTT":                                                                 //
  //=========================================================================//
  TT ToTT(TDB a_tdb)
  {
    // Get the data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    Time  diff[1];   // 1 CoOrd only
    Bits::GetCoOrds<Bits::Object::TT_TDB>(rec, dt, diff, nullptr);

    // A hack to set the TT rep directly:
    return TT() + (a_tdb.GetTimeSinceEpoch() + diff[0]);
  }

  //=========================================================================//
  // "ToTDB":                                                                //
  //=========================================================================//
  TDB ToTDB(TT a_tt)
  {
    // Solve iteratively:     tdb = tt - TT_TDB(tdb),
    // with the initial cond: tdb = tt ;
    // iterations should converge fairly quickly:
    //
    // A hack to set the TDB rep directly:
    TDB tdb(TDB() + a_tt.GetTimeSinceEpoch());

    constexpr int MaxIters = 10;
    int    i = 0;
    for (; i < MaxIters; ++i)
    {
      // Get the data for a given TDB instant. XXX: In general, if we are near
      // a node, the "rec" may change during iterations:
      auto [rec, dt] = Bits::GetRecord(tdb);

      Time  diff[1];   // 1 CoOrd only
      Bits::GetCoOrds<Bits::Object::TT_TDB>(rec, dt, diff, nullptr);

      // A hack again:
      TT    tt1(TT() + tdb.GetTimeSinceEpoch() + diff[0]);

      // For extra safety, we require a nanosecond precision:
      if (Abs(tt1 - a_tt) < Time(1e-9))
        // This "tdb" is good enough:
        break;

      // Otherwise, update "tdb" and continue. A hack again:
      tdb  = TDB() + (a_tt.GetTimeSinceEpoch() - diff[0]);
    }
    // Check for (unlikely) divergence:
    assert(i < MaxIters);
    return tdb;
  }

  //=========================================================================//
  // Earth Nutations (Long-Period Only):                                     //
  //=========================================================================//
  //                                      [d(psi),  d(eps)]
  void GetEarthNutations(TDB a_tdb, Angle a_nuts[2])
  {
    assert(a_nuts != nullptr);
    auto [rec, dt] = Bits::GetRecord(a_tdb);
    Bits::GetCoOrds <Bits::Object::EarthNutations>(rec, dt, a_nuts, nullptr);
  }

  //=========================================================================//
  // Moon Librations:                                                        //
  //=========================================================================//
  //                                      [phi, theta, psi]  ...
  void GetMoonLibrations(TDB a_tdb, Angle a_librs[3], AngVel a_libr_dots[3])
  {
    assert(a_librs != nullptr);     // "a_libr_dots" may be NULL
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    Bits::GetCoOrds <Bits::Object::MoonLibrations>
      (rec, dt, a_librs, a_libr_dots);
  }
}
// End namespace SpaceBallistics
