// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/Vector3D.hpp":                 //
//      Mechanical 3D Vectors, Parameterised by the CoOrd System (COS):      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
#include "SpaceBallistics/CoOrds/TimeScales.h"
#include <type_traits>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Vector3D" Class:                                                       //
  //=========================================================================//
  // The Vectors are  just arrays of size 3 of "DQ"s  which are assumed to be
  // instances of "DimQ"  (XXX: this is not enforced;  vectors of other types
  // are also allowed w/o any safety risks).
  // IN ADDITION, the "Vector3D" type is parameterised by the CoOrdinate System
  // ("COS") and OPTIONALLY, by the Body this Vector describes (or UNKNOWN if
  // eg it is a SpaceCraft rather than a well-known Body).
  // IMPORTANT: All Non-Inertial (Non-BaryCentric) COSes are actually SnapShots
  // taken as some Time Instant. Unfortunately, it is not possible to install
  // that Time Instant statically into the COS type,  because it is typically
  // known at run-time only. Instead, we install it the "Vector3D" as "m_cosTS".
  // In addition, "m_vecTS" is the TimeStamp of this "Vector3D" itself, is the
  // time for which this Pos, Vel etc of "B" is taken:
  //
  template<typename DQ, typename COS, Body B = Body::UNDEFINED>
  class Vector3D
  {
  public:
    //-----------------------------------------------------------------------//
    // TimeScale Used:                                                       //
    //-----------------------------------------------------------------------//
    using TS = typename COS::TimeScale;

  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    typedef DQ DQArr3[3];
    TS      m_cosTS; // COS SnapShot TimeStamp
    TS      m_vecTS; // TimeStamp of this Vector
    DQArr3  m_arr;

  public:
    //-----------------------------------------------------------------------//
    // Ctors, Assignment, Equality:                                          //
    //-----------------------------------------------------------------------//
    // Default Ctor: TimeStamp is set to UnDef, the actual Array to 0s:
    constexpr Vector3D()
    : m_cosTS(TS::UnDef()),
      m_vecTS(TS::UnDef()),
      m_arr  ()
    {}

    // Copy Ctor and assignment are auto-generated:
    constexpr Vector3D            (Vector3D const&) = default;
    constexpr Vector3D& operator= (Vector3D const&) = default;

    // Non-Default Ctors:
    constexpr Vector3D(TS a_cos_ts, TS a_vec_ts, DQ a_x, DQ a_y, DQ a_z)
    : m_cosTS(a_cos_ts),
      m_vecTS(a_vec_ts),
      m_arr { a_x, a_y, a_z }
    {}

    constexpr Vector3D(TS a_cos_ts, TS a_vec_ts, DQArr3 const& a_arr)
    : m_cosTS(a_cos_ts),
      m_vecTS(a_vec_ts),
      m_arr{  a_arr[0], a_arr[1], a_arr[2] }
    {}

    // "Init":
    constexpr void Init(TS a_cos_ts, TS a_vec_ts, DQ const a_arr[3])
    {
      m_cosTS  = a_cos_ts;
      m_vecTS  = a_vec_ts;
      m_arr[0] = a_arr[0];
      m_arr[1] = a_arr[1];
      m_arr[2] = a_arr[2];
    }

    // NB: Equality is undefined (not just "false") if the COS or Vec TimeStamps
    // do not match:
    constexpr bool operator== (Vector3D const& a_right) const
    {
      CheckTSs(a_right);
      return m_arr[0] == a_right.m_arr[0] &&
             m_arr[1] == a_right.m_arr[1] &&
             m_arr[2] == a_right.m_arr[2];
    }

    //-----------------------------------------------------------------------//
    // Accessors for Entries:                                                //
    //-----------------------------------------------------------------------//
    // NB: Entries are MUTABLE for outsize callERs (but not the TimeStamp):
    //
    // By Index:
    constexpr DQ& operator[] (int a_i)
    {
      assert(0  <= a_i && a_i < 3);
      return m_arr[a_i];
    }

    constexpr DQ const& operator[] (int a_i) const
    {
      assert(0  <= a_i && a_i < 3);
      return m_arr[a_i];
    }

    // By Name:
    constexpr DQ&        x()       { return m_arr[0]; }
    constexpr DQ const&  x() const { return m_arr[0]; }

    constexpr DQ&        y()       { return m_arr[1]; }
    constexpr DQ const&  y() const { return m_arr[1]; }

    constexpr DQ&        z()       { return m_arr[2]; }
    constexpr DQ const&  z() const { return m_arr[2]; }

    // For the TimeStamps:
    constexpr TS  GetCOSTS() const { return m_cosTS;  }
    constexpr TS  GetVecTS() const { return m_vecTS;  }

    // XXX: In some very rare cases (eg in "MechElement"), we may need writable
    // access to the COS TimeStamp. USE IT WITH EXTREME CARE!
    constexpr TS& GetCOSTS()       { return m_cosTS;  }

    //-----------------------------------------------------------------------//
    // Vector Space Operations:                                              //
    //-----------------------------------------------------------------------//
    // First of all, Unary "+" and "-":
    //
    constexpr Vector3D operator+ () const { return *this;  }

    constexpr Vector3D operator- () const
      { return  Vector3D { m_cosTS, m_vecTS, - x(), - y(), - z() }; }

    // Now the Binary ops:
    //
    // XXX:  In "+", if one of the Arg Vectors is with  UNDEFINED Body, the type
    // check is relaxed, as the result will also be for UNDEFINED.  But the COS
    // TimeStamps must match:
    //
    template<Body R>
    constexpr  Vector3D<DQ, COS, UnifyBodies(B,R)>
    operator+ (Vector3D<DQ, COS, R> const& a_right) const
    {
      return   Vector3D<DQ, COS, UnifyBodies(B,R)>
      (
        UnifyCOSTSs(a_right),
        UnifyVecTSs(a_right),
        x() + a_right.x(),
        y() + a_right.y(),
        z() + a_right.z()
      );
    }

    constexpr Vector3D& operator+= (Vector3D const& a_right)
    {
      // XXX: No unification of "Body"es is attempted here -- they must match
      // to begin with, because the "Body" ofthe LHS cannot be changed;   but
      // TimeStamps are still unified:
      //
      m_cosTS   = UnifyCOSTSs(a_right);
      m_vecTS   = UnifyVecTSs(a_right);
      m_arr[0] += a_right.m_arr[0];
      m_arr[1] += a_right.m_arr[1];
      m_arr[2] += a_right.m_arr[2];
      return *this;
    }

    // XXX:  In "-", if one of the Arg Vectors is with  UNDEFINED Body, the type
    // check is relaxed, as the result will also be for UNDEFINED.   Again, the
    // COS TimeStamps must match:
    //
    template<Body R>
    constexpr  Vector3D<DQ, COS, UnifyBodies(B,R)>
    operator- (Vector3D<DQ, COS, R> const& a_right) const
    {
      return   Vector3D<DQ, COS, UnifyBodies(B,R)>
      (
        UnifyCOSTSs(a_right),
        UnifyVecTSs(a_right),
        x() - a_right.x(),
        y() - a_right.y(),
        z() - a_right.z()
      );
    }

    constexpr Vector3D& operator-= (Vector3D const& a_right)
    {
      // XXX: No unification of "Body"es is attempted here -- they must match
      // to begin with, because the "Body" ofthe LHS cannot be changed;   but
      // TimeStamps are still unified:
      //
      m_cosTS   = UnifyCOSTSs(a_right);
      m_vecTS   = UnifyVecTSs(a_right);
      m_arr[0] -= a_right.m_arr[0];
      m_arr[1] -= a_right.m_arr[1];
      m_arr[2] -= a_right.m_arr[2];
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Multiplication by any scalar:                                         //
    //-----------------------------------------------------------------------//
    template<typename DR>
    constexpr Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, B> operator* (DR a_k)
    const
    {
      return  Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, B>
        {m_cosTS,  m_vecTS,  m_arr[0] * a_k,  m_arr[1] * a_k,  m_arr[2] * a_k};
    }

    // In-place multiplication is only possible for a "double" scalar:
    constexpr Vector3D& operator*= (double a_k)
    {
      m_arr[0] *= a_k;
      m_arr[1] *= a_k;
      m_arr[2] *= a_k;
      return *this;
    }

    template<typename DR>
    constexpr friend  Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, B> operator*
      (DR a_k, Vector3D const& a_right)
      { return a_right * a_k; }

    //-----------------------------------------------------------------------//
    // Division by any scalar:                                               //
    //-----------------------------------------------------------------------//
    template<typename DR>
    constexpr Vector3D<decltype(DQ(1.0) / DR(1.0)), COS, B> operator/ (DR a_k)
    const
    {
      assert(!::SpaceBallistics::IsZero(a_k));
      return Vector3D<decltype(DQ(1.0) / DR(1.0)), COS, B>
             {m_cosTS, m_vecTS, m_arr[0] / a_k, m_arr[1] / a_k, m_arr[2] / a_k};
    }

    // In-place division is only possible for a "double" scalar:
    constexpr Vector3D& operator/= (double a_k)
    {
      assert(a_k != 0.0);
      m_arr[0] /= a_k;
      m_arr[1] /= a_k;
      m_arr[2] /= a_k;
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Vector Length:                                                        //
    //-----------------------------------------------------------------------//
    // Implemented simply as conversion to the underlying type, as well as the
    // "EuclidNorm":
    //
    constexpr operator DQ()   const
      { return SqRt(Sqr(x()) + Sqr(y()) + Sqr(z())); }

    constexpr DQ EuclidNorm() const { return operator DQ(); }

    // "IsZero": Is it a zero-vector?
    constexpr bool IsZero()   const
    {
      return ::SpaceBallistics::IsZero(x()) &&
             ::SpaceBallistics::IsZero(y()) &&
             ::SpaceBallistics::IsZero(z());
    }

    //-----------------------------------------------------------------------//
    // Dot Product:                                                          //
    //-----------------------------------------------------------------------//
    template<typename  DR, Body R>
    constexpr decltype(DQ(1.0) * DR(1.0)) DotProd
      (Vector3D<DR, COS, R> const& a_right) const
    {
      // Check that the Bodies and TimeStamps are compatible:
      (void) UnifyBodies(B, R);
      CheckTSs(a_right);
      return x() * a_right.x() + y() * a_right.y() + z() * a_right.z();
    }

    template<typename DR, Body R>
    constexpr friend decltype(DQ(1.0) * DR(1.0)) DotProd
    (
      Vector3D<DQ, COS, B> const& a_left,
      Vector3D<DR, COS, R> const& a_right
    )
    { return a_left.DotProd(a_right); }

    //-----------------------------------------------------------------------//
    // Cross Product:                                                        //
    //-----------------------------------------------------------------------//
    template<typename DR, Body R>
    constexpr Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, UnifyBodies(B,R)>
    CrossProd(Vector3D<DR, COS, R> const& a_right) const
    {
      return
        Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, UnifyBodies(B,R)>
        {
          UnifyCOSTSs(a_right),
          UnifyVecTSs(a_right),
          y() * a_right.z() - z() * a_right.y(),
          z() * a_right.x() - x() * a_right.z(),
          x() * a_right.y() - y() * a_right.x()
        };
    }

    template<typename DR, Body R>
    constexpr friend Vector3D<decltype(DQ(1.0) * DR(1.0)), COS, B> CrossProd
    (
      Vector3D<DQ, COS, B> const& a_left,
      Vector3D<DR, COS, R> const& a_right
    )
    { return a_left.CrossProd(a_right); }

    //-----------------------------------------------------------------------//
    // Approximate Equality of Vectors:                                      //
    //-----------------------------------------------------------------------//
    template<typename DR, Body R>
    constexpr bool ApproxEquals
    (
      Vector3D<DR, COS, R> const& a_right,
      double                      a_tol = DefaultTol<double>
    )
    const
    {
      // XXX: Approximate Equality is undefined (not just "false") if COS Time-
      // Stamps do not match:
      (void) UnifyBodies(B, R);
      CheckTSs(a_right);
      return x().ApproxEquals(a_right.x(), a_tol) &&
             y().ApproxEquals(a_right.y(), a_tol) &&
             z().ApproxEquals(a_right.z(), a_tol);
    }

    //-----------------------------------------------------------------------//
    // Direct Access to the Underlying Array:                                //
    //-----------------------------------------------------------------------//
    // XXX: FOR OPTIMISATION ONLY. USE WITH CARE. KNOW WHAT YOU ARE DOING:
    //
    DQArr3 const& GetArr() const { return m_arr; }
    DQArr3&       GetArr()       { return m_arr; }

    //-----------------------------------------------------------------------//
    // TimeStamps Mgmt:                                                      //
    //-----------------------------------------------------------------------//
    // NB: "a_right" must be a Vector3D with the same COS and a compatible Body,
    // but possibly a different DimType ("DR"):
    //
    template<typename DR, Body R>
    constexpr void CheckTSs(Vector3D<DR, COS, R> const& DEBUG_ONLY(a_right))
    const
    {
      (void) UnifyBodies(B, R);
      DEBUG_ONLY
      (
        TS rCOSTS = a_right.GetCOSTS();
        TS rVecTS = a_right.GetVecTS();
      )
      assert
       ((m_cosTS.IsUnDef() || rCOSTS.IsUnDef() || m_cosTS == rCOSTS) &&
        (m_vecTS.IsUnDef() || rVecTS.IsUnDef() || m_vecTS == rVecTS));
    }

    template<typename DR, Body R>
    constexpr TS UnifyCOSTSs(Vector3D<DR, COS, R> const& a_right) const
      { return   UnifyTSs   (m_cosTS, a_right.GetCOSTS()); }

    template<typename DR, Body R>
    constexpr TS UnifyVecTSs(Vector3D<DR, COS, R> const& a_right) const
      { return   UnifyTSs   (m_vecTS, a_right.GetVecTS()); }

    //-----------------------------------------------------------------------//
    // "To_Len[_m]", "To_Len_km" Conversions:                                //
    //-----------------------------------------------------------------------//
    // Extending similar functions acting on Scalar "DimQs" to "Vector3D", just
    // for convenience:
    //
    constexpr Vector3D<decltype(::SpaceBallistics::To_Len(DQ(1.0))), COS, B>
    To_Len()  const
    {
      return
        Vector3D<decltype(::SpaceBallistics::To_Len(DQ(1.0))), COS, B>
        {
          m_cosTS,
          m_vecTS,
          ::SpaceBallistics::To_Len(x()),
          ::SpaceBallistics::To_Len(y()),
          ::SpaceBallistics::To_Len(z())
        };
    }

    constexpr Vector3D<decltype(::SpaceBallistics::To_Len(DQ(1.0))), COS, B>
    To_Len_m() const
      { return this->To_Len(); }

    constexpr Vector3D<decltype(::SpaceBallistics::To_Len_km(DQ(1.0))), COS, B>
    To_Len_km() const
    {
      return
        Vector3D<decltype(::SpaceBallistics::To_Len_km(DQ(1.0))), COS, B>
        {
          m_cosTS,
          m_vecTS,
          ::SpaceBallistics::To_Len_km(x()),
          ::SpaceBallistics::To_Len_km(y()),
          ::SpaceBallistics::To_Len_km(z())
        };
    }

    //=======================================================================//
    // "ToUnDefBody":                                                        //
    //=======================================================================//
    // Any Ref or Ptr to a Body-Specific Vector can be cast into one with
    // Body::UNDEFINED. Since "reinterpret_cast" is used, these functions
    // are NOT "constexpr":
    //
    Vector3D<DQ, COS, Body::UNDEFINED> const&   ToUnDefBody() const
    { return
        reinterpret_cast<Vector3D<DQ, COS, Body::UNDEFINED> const&> (*this); }

    Vector3D<DQ, COS, Body::UNDEFINED>&         ToUnDefBody()
      { return reinterpret_cast<Vector3D<DQ, COS, Body::UNDEFINED>&>(*this); }
  };

  //=========================================================================//
  // "ToSpecBody":                                                           //
  //=========================================================================//
  // The inverse of "ToUnDefBody": Converting a Ref or Ptr to the Vector with
  // an UnDef Body, into a Body-Specific one.  THESE FUNCTIONS SHOULD BE USED
  // WITH EXTREME CARE! They are currently only used in DE440T, for optimisa-
  // tion:
  //
  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B> const&
    ToSpecBody(Vector3D<DQ, COS, Body::UNDEFINED>  const&  a_vec)
    { return reinterpret_cast<Vector3D<DQ, COS, B> const&>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B> const*
    ToSpecBody(Vector3D<DQ, COS, Body::UNDEFINED>  const*  a_vec)
    { return reinterpret_cast<Vector3D<DQ, COS, B> const*>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B>&
    ToSpecBody(Vector3D<DQ, COS, Body::UNDEFINED>&    a_vec)
    { return reinterpret_cast<Vector3D<DQ, COS, B> &>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B>*
    ToSpecBody(Vector3D<DQ, COS, Body::UNDEFINED>*    a_vec)
    { return reinterpret_cast<Vector3D<DQ, COS, B> *>(a_vec); }

  //=========================================================================//
  // Common Mechanical "Vector3D"s:                                          //
  //=========================================================================//
  // Macro for declaring a DimQ Vector (or a diagonal Tensor). IMPORTANT: The
  // Vectors are parameterised by the CoOrd System ("COS"):
# ifdef DCL_VEC
# undef DCL_VEC
# endif
# define DCL_VEC(T) \
  template<typename COS,   Body B = Body::UNDEFINED> \
  using T##V = Vector3D<T, COS, B>;

  DCL_VEC(DimLess)  // DimLess  (eg directional vector)
  DCL_VEC(Len)      // Position Vector ("Radius-Vector")
  DCL_VEC(LenK)     // Position Vector (AstroDynamical, km)
  DCL_VEC(Vel)      // Velocity Vector
  DCL_VEC(VelK)     // Velocity Vector (AstroDynamical, km/sec)
  DCL_VEC(Acc)      // Acceleration Vector
  DCL_VEC(Force)    // Force Vector
  DCL_VEC(AngVel)   // Angular Velocity Vector
  DCL_VEC(AngAcc)   // Angular Acceleration Vector
  DCL_VEC(AngMom)   // Angular ("Kinetic") Momentum Vector
  DCL_VEC(Torq)     // Rotational Moment of Force (Torque) Vector

  // Alias for the Position vectors: "Len[K]V" -> "Pos[K]V":
  template<typename   COS, Body B = Body::UNDEFINED>
  using PosV  = LenV <COS, B>;

  template<typename   COS, Body B = Body::UNDEFINED>
  using PosKV = LenKV<COS, B>;

  // NB:  The MoI and its Rate of Change are in general not Vectors, but rather,
  // 3*3 Tensors. XXX: For the moment, we only consider those tensors in their
  // principal axes, so they have a diagonal form  and represented by Vectors:
  //
  DCL_VEC(MoI)      // Moments of Inertia
  DCL_VEC(MoIRate)  // MoI Change Rates
# undef DCL_VEC
}
// End namespace SpaceBallistics
