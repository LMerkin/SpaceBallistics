// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/Maths/RotationMatrices.hpp":                 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/Utils.hpp"

namespace SpaceBallistics
{
  //=========================================================================//
  // "Mtx33": 3*3 Matrix Type:                                               //
  //=========================================================================//
  class Mtx33
  {
  private:
    // Data:
    double m_M[3][3];

  public:
    // Default Ctor: Initialise the Mtx to 0:
    constexpr Mtx33()
    {
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        m_M[i][j] = 0.0;
    }

    // Copy Ctor, comparison and assignment are auto-generated:
    constexpr Mtx33             (Mtx33 const&)       = default;
    constexpr Mtx33& operator=  (Mtx33 const&)       = default;
    constexpr bool   operator== (Mtx33 const&) const = default;

    // Direct Access to Entries:
    constexpr double& operator() (int i, int j)
    {
      assert(0 <= i && i < 3 && 0 <= j && j < 3);
      return m_M[i][j];
    }

    constexpr double const& operator() (int i, int j) const
    {
      assert(0 <= i && i < 3 && 0 <= j && j < 3);
      return m_M[i][j];
    }
  };

  //=========================================================================//
  // "MkMtsR{1,2,3}":                                                        //
  //=========================================================================//
  // Constructs the 3D Rotation Matrices R1, R2, R3  (the naming convention is
  // taken from Astrometry), where rotation is performed for the angle "theta"
  // (positive direction is counter-clock-wise as usual) around the X, Y and Z
  // axes, resp.
  // Each function fills in 2 matrices: for the Inverse ("New-of-Old") and the
  // Direct ("Old-of-New")  transforms:
  //
  constexpr void MkMtsR1(Angle a_theta, Mtx33* a_R1, Mtx33* a_invR1)
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R1  != nullptr);
    Mtx33& R1 = *a_R1;
    R1(0,0) = 1.0;        R1(0,1) =  0.0;       R1(0,2) =  0.0;
    R1(1,0) = 0.0;        R1(1,1) =  cosTh;     R1(1,2) =  sinTh;
    R1(2,0) = 0.0;        R1(2,1) = -sinTh;     R1(2,2) =  cosTh;

    if (a_invR1 != nullptr)
    {
      Mtx33& invR1 = *a_invR1;
      invR1(0,0) = 1.0;   invR1(0,1) = 0.0;     invR1(0,2) =  0.0;
      invR1(1,0) = 0.0;   invR1(1,1) = cosTh;   invR1(1,2) = -sinTh;
      invR1(2,0) = 0.0;   invR1(2,1) = sinTh;   invR1(2,2) =  cosTh;
    }
  }

  constexpr void MkMtsR2(Angle a_theta, Mtx33* a_R2, Mtx33* a_invR2)
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R2  != nullptr);
    Mtx33& R2 = *a_R2;
    R2(0,0) = cosTh;      R2(0,1) =  0.0;       R2(0,2) = sinTh;
    R2(1,0) = 0.0;        R2(1,1) =  1.0;       R2(1,2) = 0.0;
    R2(2,0) = -sinTh;     R2(2,1) =  0.0;       R2(2,2) = cosTh;

    if (a_invR2 != nullptr)
    {
      Mtx33& invR2 = *a_invR2;
      invR2(0,0) = cosTh; invR2(0,1) =  0.0;    invR2(0,2) = -sinTh;
      invR2(1,0) = 0.0;   invR2(1,1) =  1.0;    invR2(1,2) = 0.0;
      invR2(2,0) = sinTh; invR2(2,1) =  0.0;    invR2(2,2) = cosTh;
    }
  }

  constexpr void MkMtsR3(Angle a_theta, Mtx33* a_R3, Mtx33* a_invR3)
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R3  != nullptr);
    Mtx33& R3 = *a_R3;

    R3(0,0) = cosTh;      R3(0,1) = sinTh;      R3(0,2) = 0.0;
    R3(1,0) = -sinTh;     R3(1,1) = cosTh;      R3(1,2) = 0.0;
    R3(2,0) = 0.0;        R3(2,1) =  0.0;       R3(2,2) = 1.0;

    if (a_invR3 != nullptr)
    {
      Mtx33& invR3 = *a_invR3;
      invR3(0,0) = cosTh; invR3(0,1) = -sinTh;  invR3(0,2) = 0.0;
      invR3(1,0) = sinTh; invR3(1,1) =  cosTh;  invR3(1,2) = 0.0;
      invR3(2,0) = 0.0;   invR3(2,1) =  0.0;    invR3(2,2) = 1.0;
    }
  }

  //=========================================================================//
  // 3x3 Matrix Multiplication:                                              //
  //=========================================================================//
  // XXX: Not specifically optimised; but hopefully, the compiler will optimise
  // it well enough:
  //
  constexpr void MtxMult33
  (
    Mtx33 const& a_left,
    Mtx33 const& a_right,
    Mtx33*       a_res
  )
  {
    assert(a_res != nullptr);
    Mtx33&   res  = *a_res;

    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      res(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
        res(i,j) += a_left(i,k) * a_right(k,j);
    }
  }

  //=========================================================================//
  // 3x3 Matrix-Vector Multiplication (Low-Level):                           //
  //=========================================================================//
  // NB: The Matrix is assumed to be of "double"s, but the Vector can be of any
  // suitable type (presumably a "DimQ"):
  //
  template<typename DQ>
  constexpr void MVMult33
  (
    Mtx33 const& a_mtx,
    DQ    const  a_vec[3],
    DQ           a_res[3]
  )
  {
    assert(a_vec != nullptr && a_res != nullptr);
    for (int i = 0; i < 3; ++i)
    {
      a_res[i] = DQ(0.0);
      for (int j = 0; j < 3; ++j)
        a_res[i] += a_mtx(i,j) * a_vec[j];
    }
  }
}
