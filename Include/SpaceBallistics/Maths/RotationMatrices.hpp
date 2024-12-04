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
  // "MkMtsR{1,2,3}":                                                        //
  //=========================================================================//
  // Constructs the 3D Rotation Matrices R1, R2, R3  (the naming convention is
  // taken from Astrometry), where rotation is performed for the angle "theta"
  // (positive direction is counter-clock-wise as usual) around the X, Y and Z
  // axes, resp.
  // Each function fills in 2 matrices: for the Inverse ("New-of-Old") and the
  // Direct ("Old-of-New")  transforms:
  //
  constexpr void MkMtsR1(Angle a_theta, double a_R1[3][3], double a_invR1[3][3])
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R1  != nullptr);
    a_R1   [0][0] = 1.0;     a_R1   [0][1] =  0.0;    a_R1   [0][2] =  0.0;
    a_R1   [1][0] = 0.0;     a_R1   [1][1] =  cosTh;  a_R1   [1][2] =  sinTh;
    a_R1   [2][0] = 0.0;     a_R1   [2][1] = -sinTh;  a_R1   [2][2] =  cosTh;

    if (a_invR1 != nullptr)
    {
      a_invR1[0][0] = 1.0;   a_invR1[0][1] =  0.0;    a_invR1[0][2] =  0.0;
      a_invR1[1][0] = 0.0;   a_invR1[1][1] =  cosTh;  a_invR1[1][2] = -sinTh;
      a_invR1[2][0] = 0.0;   a_invR1[2][1] =  sinTh;  a_invR1[2][2] =  cosTh;
    }
  }

  constexpr void MkMtsR2(Angle a_theta, double a_R2[3][3], double a_invR2[3][3])
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R2  != nullptr);
    a_R2   [0][0] = cosTh;   a_R2   [0][1] =  0.0;    a_R2   [0][2] = sinTh;
    a_R2   [1][0] = 0.0;     a_R2   [1][1] =  1.0;    a_R2   [1][2] = 0.0;
    a_R2   [2][0] = -sinTh;  a_R2   [2][1] =  0.0;    a_R2   [2][2] = cosTh;

    if (a_invR2 != nullptr)
    {
      a_invR2[0][0] = cosTh; a_invR2[0][1] =  0.0;    a_invR2[0][2] = -sinTh;
      a_invR2[1][0] = 0.0;   a_invR2[1][1] =  1.0;    a_invR2[1][2] = 0.0;
      a_invR2[2][0] = sinTh; a_invR2[2][1] =  0.0;    a_invR2[2][2] = cosTh;
    }
  }

  constexpr void MkMtsR3(Angle a_theta, double a_R3[3][3], double a_invR3[3][3])
  {
    double cosTh  = Cos(double(a_theta));
    double sinTh  = Sin(double(a_theta));

    assert(a_R3  != nullptr);
    a_R3   [0][0] = cosTh;   a_R3   [0][1] = sinTh;   a_R3   [0][2] = 0.0;
    a_R3   [1][0] = -sinTh;  a_R3   [1][1] = cosTh;   a_R3   [1][2] = 0.0;
    a_R3   [2][0] = 0.0;     a_R3   [2][1] =  0.0;    a_R3   [2][2] = 1.0;

    if (a_invR3 != nullptr)
    {
      a_invR3[0][0] = cosTh; a_invR3[0][1] = -sinTh;  a_invR3[0][2] = 0.0;
      a_invR3[1][0] = sinTh; a_invR3[1][1] =  cosTh;  a_invR3[1][2] = 0.0;
      a_invR3[2][0] = 0.0;   a_invR3[2][1] =  0.0;    a_invR3[2][2] = 1.0;
    }
  }

  //=========================================================================//
  // 3x3 Matrix Multiplication:                                              //
  //=========================================================================//
  // XXX: Not specifically optimised; but hopefully, the compiler will optimise
  // it well enough:
  //
  constexpr void MtxMult3
  (
    double const a_left [3][3],
    double const a_right[3][3],
    double       a_res  [3][3]
  )
  {
    assert(a_left != nullptr && a_right != nullptr && a_res != nullptr);
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      a_res[i][j] = 0.0;
      for (int k = 0; k < 3; ++k)
        a_res[i][j] += a_left[i][k] * a_right[k][j];
    }
  }

  //=========================================================================//
  // 3x3 Matrix-Vector Multiplication:                                       //
  //=========================================================================//
  // NB: The Matrix is assumed to be of "double"s, but the Vector can be of any
  // suitable type (presumably a "DimQ"):
  //
  template<typename DQ>
  constexpr void MVMult3
  (
    double const a_mtx[3][3],
    DQ     const a_vec[3],
    DQ           a_res[3]
  )
  {
    assert(a_mtx != nullptr && a_vec != nullptr && a_res != nullptr);
    for (int i = 0; i < 3; ++i)
    {
      a_res[i] = DQ(0.0);
      for (int j = 0; j < 3; ++j)
        a_res[i] += a_mtx[i][j] * a_vec[j];
    }
  }
}
