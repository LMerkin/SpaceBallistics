// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/Vector3D.hpp":                 //
//      Mechanical 3D Vectors, Parameterised by the CoOrd System (COS):      //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/Bodies.h"
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
  // eg it is a SpaceCraft rather than a well-known Body):
  //
  template<typename DQ, typename COS, Body B = Body::UNDEFINED>
  class Vector3D
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    typedef DQ DQArr3[3];
    DQArr3  m_arr;   // Initialised to 0s by default

  public:
    //-----------------------------------------------------------------------//
    // Ctors, Assignment, Equality:                                          //
    //-----------------------------------------------------------------------//
    // Default Ctor, Copy Ctor, Assignment and Equality:are auto-generated; us-
    // ing the Default Ctor of "DQ":
    //
    constexpr Vector3D            ()                      = default;
    constexpr Vector3D            (Vector3D const&)       = default;
    constexpr Vector3D& operator= (Vector3D const&)       = default;
    constexpr bool      operator==(Vector3D const&) const = default;

    // Non-Default Ctors:
    constexpr Vector3D(DQ a_x, DQ a_y, DQ a_z)
    : m_arr{ a_x, a_y, a_z }
    {}

    constexpr Vector3D(DQArr3 const& a_arr)
    : m_arr{ a_arr[0], a_arr[1], a_arr[2] }
    {}

    //-----------------------------------------------------------------------//
    // Accessors for Entries:                                                //
    //-----------------------------------------------------------------------//
    // NB: Entries are in general MUTABLE:
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
    constexpr DQ&       x()       { return m_arr[0]; }
    constexpr DQ const& x() const { return m_arr[0]; }

    constexpr DQ&       y()       { return m_arr[1]; }
    constexpr DQ const& y() const { return m_arr[1]; }

    constexpr DQ&       z()       { return m_arr[2]; }
    constexpr DQ const& z() const { return m_arr[2]; }

    //-----------------------------------------------------------------------//
    // Vector Space Operations:                                              //
    //-----------------------------------------------------------------------//
    // XXX:  In "+", if one of the Arg Vectors is with  UNDEFINED Body, the type
    // check is relaxed, as the result will also be for UNDEFINED:
    //
    template<Body R>
    constexpr  Vector3D<DQ, COS, UnifyBodies(B,R)>
    operator+ (Vector3D<DQ, COS, R> const& a_right) const
    {
      return   Vector3D<DQ, COS, UnifyBodies(B,R)>
      (
        x() + a_right.x(),
        y() + a_right.y(),
        z() + a_right.z()
      );
    }

    constexpr Vector3D& operator+= (Vector3D const& a_right)
    {
      // But here no unification is attempted:
      m_arr[0] += a_right.m_arr[0];
      m_arr[1] += a_right.m_arr[1];
      m_arr[2] += a_right.m_arr[2];
      return *this;
    }

    // XXX:  In "-", if one of the Arg Vectors is with  UNDEFINED Body, the type
    // check is relaxed, as the result will also be for UNDEFINED:
    //
    template<Body R>
    constexpr  Vector3D<DQ, COS, UnifyBodies(B,R)>
    operator- (Vector3D<DQ, COS, R> const& a_right) const
    {
      return   Vector3D<DQ, COS, UnifyBodies(B,R)>
      (
        x() - a_right.x(),
        y() - a_right.y(),
        z() - a_right.z()
      );
    }

    constexpr Vector3D& operator-= (Vector3D const& a_right)
    {
      // But here no unification is attempted:
      m_arr[0] -= a_right.m_arr[0];
      m_arr[1] -= a_right.m_arr[1];
      m_arr[2] -= a_right.m_arr[2];
      return *this;
    }

    //-----------------------------------------------------------------------//
    // Multiplication by any scalar:                                         //
    //-----------------------------------------------------------------------//
    template<typename DT>
    constexpr Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B> operator* (DT a_k)
    const
    {
      return Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B>
             (m_arr[0] * a_k, m_arr[1] * a_k, m_arr[2] * a_k);
    }

    // In-place multiplication is only possible for a "double" scalar:
    constexpr Vector3D& operator*= (double a_k)
    {
      m_arr[0] *= a_k;
      m_arr[1] *= a_k;
      m_arr[2] *= a_k;
      return *this;
    }

    template<typename DT>
    constexpr friend  Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B> operator*
      (DT a_k, Vector3D const& a_right)
      { return a_right * a_k; }

    //-----------------------------------------------------------------------//
    // Division by any scalar:                                               //
    //-----------------------------------------------------------------------//
    template<typename DT>
    constexpr Vector3D<decltype(DQ(1.0) / DT(1.0)), COS, B> operator/ (DT a_k)
    const
    {
      assert(!IsZero(a_k));
      return Vector3D<decltype(DQ(1.0) / DT(1.0)), COS, B>
             (m_arr[0] / a_k, m_arr[1] / a_k, m_arr[2] / a_k);
    }

    // Overloading (NOT Specialisation) of the above:
    // Division by a "DQ" will produce a vector of "double"s:
    //
    constexpr Vector3D<double, COS, B> operator/ (DQ a_k) const
    {
      assert(!IsZero(a_k));
      return Vector3D<double, COS, B>
             (double(m_arr[0] / a_k), double(m_arr[1] / a_k),
              double(m_arr[2] / a_k));
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
    // "EuclNorm":
    //
    constexpr operator DQ() const
      { return SqRt(Sqr(x()) + Sqr(y()) + Sqr(z())); }

    constexpr DQ EuclNorm() const { return operator DQ(); }

    //-----------------------------------------------------------------------//
    // Dot Product:                                                          //
    //-----------------------------------------------------------------------//
    template<typename  DT>
    constexpr decltype(DQ(1.0) * DT(1.0)) DotProd
      (Vector3D<DT, COS, B> const& a_right) const
    { return x() * a_right.x() + y() * a_right.y() + z() * a_right.z(); }

    template<typename DT>
    constexpr friend decltype(DQ(1.0) * DT(1.0)) DotProd
    (
      Vector3D<DQ, COS, B> const& a_left,
      Vector3D<DT, COS, B> const& a_right
    )
    { return a_left.DotProd(a_right); }

    //-----------------------------------------------------------------------//
    // Cross Product:                                                        //
    //-----------------------------------------------------------------------//
    template<typename DT>
    constexpr Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B> CrossProd
      (Vector3D<DT, COS, B> const& a_right) const
    {
      return
        Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B>
        {
          y() * a_right.z() - z() * a_right.y(),
          z() * a_right.x() - x() * a_right.z(),
          x() * a_right.y() - y() * a_right.x()
        };
    }

    template<typename DT>
    constexpr friend Vector3D<decltype(DQ(1.0) * DT(1.0)), COS, B> CrossProd
    (
      Vector3D<DQ, COS, B> const& a_left,
      Vector3D<DT, COS, B> const& a_right
    )
    { return a_left.CrossProd(a_right); }

    //-----------------------------------------------------------------------//
    // Approximate Equality of Vectors:                                      //
    //-----------------------------------------------------------------------//
    constexpr bool ApproxEquals
    (
      Vector3D const& a_right,
      double          a_tol = DefaultTol<double>
    )
    const
    {
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
  };

  //=========================================================================//
  // "ToGeneric":                                                            //
  //=========================================================================//
  // Any Ref or Ptr to a Body-Specific Vector can be cast into the Generic one
  // (with Body::UNDEFINED). Since "reinterpret_casr" is used, these functions
  // are NOT "constexpr":
  //
  template<typename DQ,  typename COS, Body B>
  Vector3D<DQ, COS, Body::UNDEFINED> const&
    ToGeneric(Vector3D<DQ, COS, B>   const& a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B> const&>(a_vec); }

  template<typename DQ,  typename COS, Body B>
  Vector3D<DQ, COS, Body::UNDEFINED> const*
    ToGeneric(Vector3D<DQ, COS, B>   const* a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B> const*>(a_vec); }

  template<typename DQ,  typename COS, Body B>
  Vector3D<DQ, COS, Body::UNDEFINED>&
    ToGeneric(Vector3D<DQ, COS, B>  & a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B>&>(a_vec); }

  template<typename DQ,  typename COS, Body B>
  Vector3D<DQ, COS, Body::UNDEFINED>*
    ToGeneric(Vector3D<DQ, COS, B>  * a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B>*>(a_vec); }

  //=========================================================================//
  // "ToSpecific":                                                           //
  //=========================================================================//
  // The inverse of "ToGeneric": Converting a Ref or Ptr to the Generic Vector,
  // into a Body-Specific one. USE WITH EXTREME CARE!
  //
  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B> const&
    ToSpecific(Vector3D<DQ, COS, Body::UNDEFINED> const&  a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B>  const&>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B> const*
    ToSpecific(Vector3D<DQ, COS, Body::UNDEFINED> const*  a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B>  const*>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B>&
    ToSpecific(Vector3D<DQ, COS, Body::UNDEFINED>&  a_vec)
  { return reinterpret_cast<Vector3D<DQ, COS, B> &>(a_vec); }

  template<Body B,  typename DQ,  typename COS>
  Vector3D<DQ, COS, B>*
    ToSpecific(Vector3D<DQ, COS, Body::UNDEFINED>*  a_vec)
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
