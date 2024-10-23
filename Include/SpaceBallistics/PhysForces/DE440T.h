// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysForces/DE440T.h":                   //
//                         JPL DE440T Ephemerides                            //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <utility>

namespace SpaceBallistics::DE440T
{
  //=========================================================================//
  // "Bits": Internal Implementation:                                        //
  //=========================================================================//
  namespace Bits
  {
    //-----------------------------------------------------------------------//
    // Consts:                                                               //
    //-----------------------------------------------------------------------//
    // The Temporal Range of Records (in TDB):
    constexpr TDB From = TDB(Time_day(2'323'696.5)); // 1649-12-18.0
    constexpr TDB To   = TDB(Time_day(2'506'320.5)); // 2149-12-21.0
  
    // The Number of Records:
    constexpr int NR   = 5707;
  
    // Temporal Span of each Record:
    constexpr Time      RecSpan =  To_Time(32.0_day);
    static_assert(double(NR) * RecSpan == To - From);
  
    // Size of each record in "double"s:
    constexpr int ND   = 1122;
  
    //-----------------------------------------------------------------------//
    // Objects Provided by DE440T:                                           //
    //-----------------------------------------------------------------------//
    // Unless otherwise stated, for a given Object, DE440T provide 3D coords
    // [X, Y, Z] in the "BaryCentricCOS":
    //
    enum class Object: int
    {
      Mercury         = 0,
      Venus           = 1,
      EMB             = 2,     // Earth-Moon System BaryCenter
      Mars            = 3,
      Jupiter         = 4,
      Saturn          = 5,
      Uranus          = 6,
      Neptune         = 7,
      Pluto           = 8,
      MoonGeoC        = 9,     // Moon [X, Y, Z] in the "GeoCentricFixedCOS"
      Sun             = 10,
      EarthNutations  = 11,    // [d(psi),  d(eps)]
      MoonLibrations  = 12,    // [phi, theta, psi]
      TT_TDB          = 13     // [TT-TDB]
    };
    constexpr int NObjs = int(Object::TT_TDB) + 1;
  
    //-----------------------------------------------------------------------//
    // Chebyshev Coeffs Layout for each Object:                              //
    //-----------------------------------------------------------------------//
    // Dim1 (Outer): The Number of SubPeriods per RecSpan:
    constexpr int NSPs[NObjs]
      {  4,  2,  2,  1,  1,  1,  1,  1,  1,  8,  2,  4,  4,  8 };
  
    // Same, as a templated value:
    template<Object Obj>
    constexpr int NSP = NSPs[int(Obj)];
  
    // Dim2 (Mid):   The Number of CoOrds:
    constexpr int NCOs[NObjs]
      {  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  2,  3,  1 };
  
    // Same, as a templated value:
    template<Object Obj>
    constexpr int NCO = NCOs[int(Obj)];
  
    // Dim3 (Inner): The Number of Chebyshev Coeffs per CoOrd:
    constexpr int NCCs[NObjs]
      { 14, 10, 13, 11,  8,  7,  6,  6,  6, 13, 11, 10, 10, 13 };
  
    // Same, as a templated value:
    template<Object Obj>
    constexpr int NCC = NCCs[int(Obj)];
  
    //-----------------------------------------------------------------------//
    // Types of Chebyshev Coeffs for different "Object"s:                    //
    //-----------------------------------------------------------------------//
    // In most cases, it is "Len_km", but there are some exceptions:
    template<Object Obj>
    struct CTW                         { using T = Len_km; };

    template<>
    struct CTW<Object::EarthNutations> { using T = Angle;  };

    template<>
    struct CTW<Object::MoonLibrations> { using T = Angle;  };

    template<>
    struct CTW<Object::TT_TDB>         { using T = Time;   };
  
    // The "high-level" short-cut:
    template<Object Obj>
    using CoeffsT = CTW<Obj>::T;
  
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
      CoeffsT<Object::ObjName> const m_##ObjName \
        [NSP<Object::ObjName>]  \
        [NCO<Object::ObjName>]  \
        [NCC<Object::ObjName>];
  
      // So: Chebyshev Coeffs for all Objects:
      MK_DE440T_REC_ENTRY(Mercury)
      MK_DE440T_REC_ENTRY(Venus)
      MK_DE440T_REC_ENTRY(EMB)
      MK_DE440T_REC_ENTRY(Mars)
      MK_DE440T_REC_ENTRY(Jupiter)
      MK_DE440T_REC_ENTRY(Saturn)
      MK_DE440T_REC_ENTRY(Uranus)
      MK_DE440T_REC_ENTRY(Neptune)
      MK_DE440T_REC_ENTRY(Pluto)
      MK_DE440T_REC_ENTRY(MoonGeoC)
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
    CoeffsT <Obj> const* GetChebyshevCoeffs(Record const* a_rec, int a_s);
  
#   ifdef  MK_DE440T_GET_COEFFS
#   undef  MK_DE440T_GET_COEFFS
#   endif
#   define MK_DE440T_GET_COEFFS(ObjName)       \
    template<> \
    inline CoeffsT<Object::ObjName> const*     \
    GetChebyshevCoeffs<Object::ObjName>(Record const* a_rec, int a_s)     \
    { \
      assert(a_rec != nullptr && 0 <= a_s && a_s < NSP<Object::ObjName>); \
      return &(a_rec->m_##ObjName[a_s][0][0]); \
    }
    MK_DE440T_GET_COEFFS(Mercury)
    MK_DE440T_GET_COEFFS(Venus)
    MK_DE440T_GET_COEFFS(EMB)
    MK_DE440T_GET_COEFFS(Mars)
    MK_DE440T_GET_COEFFS(Jupiter)
    MK_DE440T_GET_COEFFS(Saturn)
    MK_DE440T_GET_COEFFS(Uranus)
    MK_DE440T_GET_COEFFS(Neptune)
    MK_DE440T_GET_COEFFS(Pluto)
    MK_DE440T_GET_COEFFS(MoonGeoC)
    MK_DE440T_GET_COEFFS(Sun)
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
    // TDB -> (Record, TimeOffset witin that Record):
    std::pair<Record const*, Time> GetRecord(TDB a_t);
  
    // (Object, Record, TimeOffset within that Record) -> CoOrds:
    template<Object Obj>
    void GetCoOrds
    (
      Record  const* a_record,
      Time           a_dt,
      CoeffsT<Obj>   (&a_res)[NCO<Obj>]   // Output for all CoOrds
    );
  }
  // End namespace "Bits"

  //=========================================================================//
  // "External" Functions:                                                   //
  //=========================================================================//
  //-------------------------------------------------------------------------//
  // TDB <-> TT Conversions:                                                 //
  //-------------------------------------------------------------------------//
  TT  TTofTDB(TDB a_tdb);
  TDB TDBofTT(TT  a_tt);

  //-----------------------------------------------------------------------//
  // "SelfTest" (for Temporal Continuity of Data):                         //
  //-----------------------------------------------------------------------//
  void SelfTest();
}
// End namespace SpaceBallistics::DE440T
