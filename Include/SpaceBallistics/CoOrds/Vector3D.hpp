// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/CoOrds/Vector3D.hpp":                 //
//      Mechanical 3D Vectors, Parameterised by the CoOrd System (COS):      //
//===========================================================================//
#pragma once

namespace SpaceBallistics
{
  //=========================================================================//
  // "Vector3D" Class:                                                       //
  //=========================================================================//
  // The Vectors are  just arrays of size 3 of "DQ"s  which are assumed to be
  // instances of "DimQ"  (XXX: this is not enforced;  vectors of other types
  // are also allowed w/o any safety risks). IN ADDITION, the "Vector3D" type
  // is parameterised by the CoOrdinate System ("COS"):
  //
  template<typename DQ, typename COS>
  class Vector3D
  {
  private:
    //-----------------------------------------------------------------------//
    // Data Flds:                                                            //
    //-----------------------------------------------------------------------//
    typedef DQ DQArr3[3];
    DQArr3  m_arr;

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

    // Non-Default Ctor:
    constexpr Vector3D(DQ a_x, DQ a_y, DQ a_z)
    : m_arr{ a_x, a_y, a_z }
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
    constexpr Vector3D  operator+  (Vector3D const& a_right) const
    {
      return  Vector3D(m_arr[0] + a_right.m_arr[0],
                       m_arr[1] + a_right.m_arr[1],
                       m_arr[2] + a_right.m_arr[2]);
    }

    constexpr Vector3D& operator+= (Vector3D const& a_right)
    {
      m_arr[0] += a_right.m_arr[0];
      m_arr[1] += a_right.m_arr[1];
      m_arr[2] += a_right.m_arr[2];
      return *this;
    }

    constexpr Vector3D  operator-  (Vector3D const& a_right) const
    {
      return  Vector3D(m_arr[0] - a_right.m_arr[0],
                       m_arr[1] - a_right.m_arr[1],
                       m_arr[2] - a_right.m_arr[2]);
    }

    constexpr Vector3D& operator-= (Vector3D const& a_right)
    {
      m_arr[0] -= a_right.m_arr[0];
      m_arr[1] -= a_right.m_arr[1];
      m_arr[2] -= a_right.m_arr[2];
      return *this;
    }

    constexpr Vector3D  operator*  (double a_k) const
      { return  Vector3D(a_k * m_arr[0], a_k * m_arr[1], a_k * m_arr[2]); }

    constexpr Vector3D& operator*= (double a_k)
    {
      m_arr[0] *= a_k;
      m_arr[1] *= a_k;
      m_arr[2] *= a_k;
      return *this;
    }

    constexpr friend Vector3D operator* (double a_k, Vector3D const& a_right)
      { return a_right * a_k; }

    // Vector Length:
    // Implemented simply as conversion to the underlying type:
    //
    constexpr operator DQ() const
      { return SqRt(Sqr(m_arr[0]) + Sqr(m_arr[1]) + Sqr(m_arr[2])); }

    //-----------------------------------------------------------------------//
    // Direct Access to the Underlying Array:                                //
    //-----------------------------------------------------------------------//
    // XXX: FOR OPTIMISATION ONLY. USE WITH CARE. KNOW WHAT YOU ARE DOING:
    //
    DQArr3 const& GetArr() const { return m_arr; }
    DQArr3&       GetArr()       { return m_arr; }
  };
}
// End namespace SpaceBallistics
