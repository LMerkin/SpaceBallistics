// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/PhysForces/DE440T.hpp":                  //
//                   Functions (but not Data) for DE440T                     //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/PhysForces/DE440T.h"
#include "SpaceBallistics/Maths/Chebyshev.hpp"
#include <cmath>
#include <cassert>

namespace SpaceBallistics::DE440T
{
  namespace Bits
  {
    //=======================================================================//
    // "GetRecord":                                                          //
    //=======================================================================//
    // Internal Util. Returns the pair:
    //   (Record Ptr corresponding to the given TDB instant,
    //    Time Offset within the Record):
    //
    std::pair<Record const*, Time> GetRecord(TDB a_tdb)
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
      constexpr int  NSP = NSPs[int(Obj)];
      constexpr Time SP  = RecSpan     / double(NSP);  // SubPeriod Length
      double         s0  = double(a_dt / SP);
      double         sf  = std::floor(s0);
      int            s   = int       (sf);
      assert( 0 <=   s   &&  s    <  NSP);

      // Compute the offset "tau" from the beginning of the SubPeriod:
      Time   tau  = a_dt - sf * SP;
      assert(tau <= a_dt);

      // In any case, we must get valid "s" and "tau":
      assert(0 <= s && s < NSP && !IsNeg(tau) && tau < SP && tau <= a_dt);

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
        a_pos[i] =   Chebyshev::Sum1T<double, CT<Obj>>(NCC<Obj>-1, coeffsI, x)
                     + 0.5 * coeffsI[0];

        if (a_vel != nullptr)
          // Provide the Dots (Velocities) as well. They are NOT stored  in the
          // DE440T data sets, so we have to differentiate the Chebyshev Series
          // (XXX:  do they really provide good approximations in C1,  not just
          // in C0???):
          a_vel[i] = Chebyshev::SumDT<double, CT<Obj>>(NCC<Obj>-1, coeffsI, x)
                   * dxdt;
      }
      // All Done!
    }
  }
  // End namespace Bits

  //=========================================================================//
  // "GetPlanet[s]BPV[s]":                                                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // With Compile-Time Object Selection:                                     //
  //-------------------------------------------------------------------------//
  template<Body BodyName>
  void GetPlanetBPV
  (
    TDB        a_tdb,
    PosKVBary* a_pos,    // Output (Position)
    VelKVBary* a_vel     // Output (Velocity); may be NULL
  )
  {
    // This is the generic case for the Sun and Planets except the Earth; the
    // Moon is also not supported here (as we seldom need its BaryCentruc PV):
    static_assert(BodyName != Body::Earth && BodyName != Body::Moon);
    assert(a_pos != nullptr);

    // Get the record data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    // Translate "Body" into "Object", they use consistent enumeration:
    constexpr Bits::Object Obj = Bits::Object(int(BodyName));

    // Store the result ditrectly in the underlying arrays of "a_pos", "a_vel":
    Bits::GetCoOrds<Obj>
    (
      rec,
      dt,
      a_pos->GetArr(),
      a_vel != nullptr ? a_vel->GetArr() : nullptr
    );
  }

  //-------------------------------------------------------------------------//
  // Specialisation for the Earth:                                           //
  //-------------------------------------------------------------------------//
  template<>
  void GetPlanetBPV<Body::Earth>
  (
    TDB        a_tdb,
    PosKVBary* a_pos,    // Output (Position)
    VelKVBary* a_vel     // Output (Velocity); may be NULL
  )
  {
    assert(a_pos != nullptr);

    // Get the record data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    // First, get the PV of the Earth-Moon System BaryCenter:
    PosKVBary posEMB;
    VelKVBary velEMB;
    Bits::GetCoOrds<Bits::Object::EMB>
    (
      rec,
      dt,
      posEMB.GetArr(),
      (a_vel != nullptr) ? velEMB.GetArr()  : nullptr
    );

    // Now get the GeoCentric PV of the Moon:
    PosKVGeoF posMoon;
    VelKVGeoF velMoon;
    Bits::GetCoOrds<Bits::Object::Moon>
    (
      rec,
      dt,
      posMoon.GetArr(),
      (a_vel != nullptr) ? velMoon.GetArr() : nullptr
    );

    // The GeoCentric  PV of the EMB: Proportional to those of the Moon:
    constexpr double mu = 1.0 / (Bits::EMRat + 1.0);

    // So  the BaryCentric PV of the Earth will be:
    (*a_pos)[0]   = posEMB[0] - mu * posMoon[0];
    (*a_pos)[1]   = posEMB[1] - mu * posMoon[1];
    (*a_pos)[2]   = posEMB[2] - mu * posMoon[2];

    if (a_vel != nullptr)
    {
      (*a_vel)[0] = velEMB[0] - mu * velMoon[0];
      (*a_vel)[1] = velEMB[1] - mu * velMoon[1];
      (*a_vel)[2] = velEMB[2] - mu * velMoon[2];
    }
  }

  //-------------------------------------------------------------------------//
  // With Run-Time Body Selection:                                           //
  //-------------------------------------------------------------------------//
  void GetPlanetBPV
  (
    Body       a_body,
    TDB        a_tdb,
    PosKVBary* a_pos,    // Output (Position)
    VelKVBary* a_vel     // Output (Velocity); may be NULL
  )
  {
    switch (a_body)
    {
#   ifdef  DE440T_BODY_CASE
#   undef  DE440T_BODY_CASE
#   endif
#   define DE440T_BODY_CASE(BodyName) \
      case Body::BodyName: \
           GetPlanetBPV<Body::BodyName>(a_tdb, a_pos, a_vel); \
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
      DE440T_BODY_CASE(Pluto)
#   undef  DE440T_BODY_CASE
    default:
      assert(false);
    }
  }

  //-------------------------------------------------------------------------//
  // For the 10 Objects Simultaneously:                                      //
  //-------------------------------------------------------------------------//
  // Since the primary use case of this function is to compute the planetary po-
  // sitions for integration of Minor Solar System Bodies' Orbits,  it produces
  // the co-ords of EMB rather than Earth and Moon separately. The output array
  // is (same enumeration as that of "Object"s):
  //
  // [0: Sun,    1: Mercury, 2: Venus,   3: EMB, 4: Mars, 5: Jupiter,
  //  6: Saturn, 7: Uranus,  8: Neptune, 9: Pluto]:
  //
  void GetPlanetsBPVs
  (
    TDB        a_tdb,
    PosKVBary  a_poss[10], // Output
    VelKVBary  a_vels[10]  // Output (again, may be NULL)
  )
  {
    // Optimisation: The Record is fetched just once for all Objects:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

#   ifdef  DE440T_OBJ_PV
#   undef  DE440T_OBJ_PV
#   endif
#   define DE440T_OBJ_PV(ObjName) \
    { \
      Bits::GetCoOrds<Bits::Object::ObjName> \
        (rec, \
         dt,  \
         a_poss  [int(Bits::Object::ObjName)].GetArr(), \
         ( a_vels != nullptr ) \
         ? a_vels[int(Bits::Object::ObjName)].GetArr()  \
         : nullptr \
      ); \
    }
    DE440T_OBJ_PV(Sun)
    DE440T_OBJ_PV(Mercury)
    DE440T_OBJ_PV(Venus)
    DE440T_OBJ_PV(EMB)
    DE440T_OBJ_PV(Mars)
    DE440T_OBJ_PV(Jupiter)
    DE440T_OBJ_PV(Saturn)
    DE440T_OBJ_PV(Uranus)
    DE440T_OBJ_PV(Neptune)
    DE440T_OBJ_PV(Pluto)
#   undef DE440T_OBJ_PV
  }

  //=========================================================================//
  // "TTofTDB":                                                              //
  //=========================================================================//
  TT TTofTDB(TDB a_tdb)
  {
    // Get the data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    Time  diff[1];   // 1 CoOrd only
    Bits::GetCoOrds<Bits::Object::TT_TDB>(rec, dt, diff);

    return TT() + (a_tdb.GetTime() + diff[0]);
  }

  //=========================================================================//
  // "TDBofTT":                                                              //
  //=========================================================================//
  TDB TDBofTT(TT a_tt)
  {
    // Solve iteratively:     tdb = tt - TT_TDB(tdb),
    // with the initial cond: tdb = tt ;
    // iterations should converge fairly quickly:
    //
    TDB tdb(TDB() + a_tt.GetTime());

    constexpr int MaxIters = 10;
    int    i = 0;
    for (; i < MaxIters; ++i)
    {
      // Get the data for a given TDB instant. XXX: In general, if we are near
      // a node, the "rec" may change during iterations:
      auto [rec, dt] = Bits::GetRecord(tdb);

      Time  diff[1];   // 1 CoOrd only
      Bits::GetCoOrds<Bits::Object::TT_TDB>(rec, dt, diff);

      TT    tt1(TT() + tdb.GetTime()  + diff[0]);

      // For extra safety, we require a nanosecond precision:
      if (Abs(tt1 - a_tt) < Time(1e-9))
        // This "tdb" is good enough:
        break;

      // Otherwise, update "tdb" and continue:
      tdb  = TDB()   + (a_tt.GetTime() - diff[0]);
    }
    // Check for (unlikely) divergence:
    assert(i < MaxIters);
    return tdb;
  }

  //=========================================================================//
  // "SelfTest":                                                             //
  //=========================================================================//
  void SelfTest()
  {
    // Verify the Temporal Continuity of DE440T Data:
    TDB expFrom = Bits::From;
    TDB expTo;    // Empty as yet

    for (int r = 0; r < Bits::NR; ++r)
    {
      Bits::Record const* rec =
        reinterpret_cast<Bits::Record const*>(&(Bits::Data[r][0]));

      TDB from(rec->m_From);
      TDB to  (rec->m_To);
      expTo   = from + Bits::RecSpan;

      assert(from == expFrom && to == expTo);
      expFrom = to;
    }
    assert(expTo == Bits::To);
  }
}
// End namespace SpaceBallistics
