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
      assert(From <= a_tdb && a_tdb <= To);

      // XXX: Be careful to prevent accidential big rounding-down errors!
      // If "a_tdb" is close enough to a node, return that node, otherwise take
      // the Floor:
      Time   off = a_tdb - From;
      assert(!IsNeg(off));

      double  r0  = double(off / RecSpan);
      double  r1  = std::round(r0);
      Time    dt;   // Initially 0
      int     r   = INT_MAX;

      if (UNLIKELY(Abs(r0 - r1) < 5e-12))
      {
        r  = int(r1);
        assert(IsZero(dt));

        // In this case, we may get r==NR, so adjust it:
        assert(0 <= r && r <= NR);
        if (r == NR)
        {
          --r;
          dt = RecSpan;
        }
      }
      else
      {
        double   r2 = std::floor(r0);
        r  = int(r2);
        dt = off - r2 * RecSpan;

        // In this case, the following must hold:
        assert(0 <= r && r < NR && IsPos(dt) && dt <  RecSpan);
      }

      // Once again: The over-all check:
      assert(0  <= r && r < NR && !IsNeg(dt));
      assert(dt < RecSpan      || (dt == RecSpan && r == NR-1));

      // We can now extract the Record:
      Record const* rec = reinterpret_cast<Record const*>(&(Data[r][0]));
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
      assert(a_record != nullptr && !IsNeg(a_dt) && a_dt <= RecSpan &&
             a_pos    != nullptr);

      // Initially "a_dt" is the temporal offset within the whole RecSpan; find
      // the corresp SubPeriod Index ("s1"):
      constexpr int  NSP = NSPs[int(Obj)];
      constexpr Time SP  = RecSpan     / double(NSP);  // SubPeriod
      double         s0  = double(a_dt / SP);
      double         s1  = std::floor(s0);
      int            s   = int       (s1);
      assert( 0 <=   s   &&  s    <= NSP);
      assert((s == NSP)  == (a_dt == RecSpan));        // Special case...

      // Compute the offset "tau" from the beginning of the SubPeriod:
      Time tau  =   a_dt - s1 * SP;

      // In theory, we must have 0 <= tau <= SP, but there may be rounding
      // errors, so enforce this constraint:
      tau = std::min(std::max(tau, 0.0_sec), SP);

      // NB: s==NSP can only happen when a_dt == RecSpan (asserted above), in
      // which case we need an adjustment:
      if (UNLIKELY(s == NSP))
      {
        // This means that "tau" is 0, but again, beware of rounding errors:
        assert(tau.ApproxEquals(0.0_sec));

        // Adjust "s" and "tau" to the previous sub-period (always exists):
        --s;
        assert(s >= 0);
        tau = SP;
      }
      // In any case, we must get valid "s" and "tau":
      assert(0 <= s && s < NSP && !IsNeg(tau) && tau <= SP);

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
        CT<Obj> const* coeffsI = coeffsBase     + i * NCC<Obj>;
        a_pos[i] =   Chebyshev::Sum1T<double, CT<Obj>>(NCC<Obj>-1, coeffsI, x);

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
  // "GetBaryPosVel":                                                        //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // With Compile-Time Object Selection:                                     //
  //-------------------------------------------------------------------------//
  template<Object Obj>
  void GetBaryPosVel
  (
    TDB        a_tdb,
    PosKVBary* a_pos,    // Output (Position)
    VelKVBary* a_vel     // Output (Velocity); may be NULL
  )
  {
    // This function is applicable to the following Objects only:
    static_assert(Obj == Object::Mercury || Obj == Object::Venus   ||
                  Obj == Object::EMB     || Obj == Object::Mars    ||
                  Obj == Object::Jupiter || Obj == Object::Saturn  ||
                  Obj == Object::Uranus  || Obj == Object::Neptune ||
                  Obj == Object::Pluto   || Obj == Object::Sun);
    assert(a_pos != nullptr);

    // Get the data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    // Store the result ditrectly in the underlying arrays of "a_pos", "a_vel":
    GetCoOrds<Obj>
      (rec, dt, a_pos->GetArr(), a_vel != nullptr ? a_vel->GetArr() : nullptr);
  }

  //-------------------------------------------------------------------------//
  // With Run-Time Object Selection:                                         //
  //-------------------------------------------------------------------------//
  void GetBaryPosVel
  (
    Object     a_obj,
    TDB        a_tdb,
    PosKVBary* a_pos,    // Output (Position)
    VelKVBary* a_vel     // Output (Velocity); may be NULL
  )
  {
    switch (a_obj)
    {
#   ifdef  DE440T_OBJ_CASE
#   undef  DE440T_OBJ_CASE
#   endif
#   define DE440T_OBJ_CASE(ObjName) \
      case Object::ObjName: \
           GetBaryPosVel<Object::ObjName>(a_tdb, a_pos, a_vel); \
           break;
      DE440T_OBJ_CASE(Mercury)
      DE440T_OBJ_CASE(Venus)
      DE440T_OBJ_CASE(EMB)
      DE440T_OBJ_CASE(Mars)
      DE440T_OBJ_CASE(Jupiter)
      DE440T_OBJ_CASE(Saturn)
      DE440T_OBJ_CASE(Uranus)
      DE440T_OBJ_CASE(Neptune)
      DE440T_OBJ_CASE(Pluto)
      DE440T_OBJ_CASE(Sun)
#   undef  DE440T_OBJ_CASE
    default:
      assert(false);
    }
  }

  //-------------------------------------------------------------------------//
  // Forr All 10 Objects Simultaneously:                                     //
  //-------------------------------------------------------------------------//

  //=========================================================================//
  // "TTofTDB":                                                              //
  //=========================================================================//
  TT TTofTDB(TDB a_tdb)
  {
    // Get the data for a given TDB instant:
    auto [rec, dt] = Bits::GetRecord(a_tdb);

    Time  diff[1];  // 1 CoOrd only
    Bits::GetCoOrds<Object::TT_TDB>(rec, dt, diff);

    return TT(a_tdb.GetTime() + diff[0]);
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
