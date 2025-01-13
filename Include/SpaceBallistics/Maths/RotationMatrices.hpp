// vim:ts=2:et
//===========================================================================//
//             "SpaceBallistics/Maths/RotationMatrices.hpp":                 //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"

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

    //-----------------------------------------------------------------------//
    // Matrix Transposition:                                                 //
    //-----------------------------------------------------------------------//
    constexpr Mtx33 Transpose() const
    {
      Mtx33 res;
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        res(i,j) = (*this)(j,i);
      return res;
    }

    //-----------------------------------------------------------------------//
    // Matrix Multiplication:                                                //
    //-----------------------------------------------------------------------//
    // XXX: Not specifically optimised; but hopefully, the compiler will optim-
    // ise  it well enough:
    //
    constexpr Mtx33 operator*(Mtx33 const& a_right) const
    {
      Mtx33 res;
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        res(i,j) = 0.0;
        for (int k = 0; k < 3; ++k)
          res(i,j) += (*this)(i,k) * a_right(k,j);
      }
      return res;
    }

    //=======================================================================//
    // "MkR[123]":                                                           //
    //=======================================================================//
    // Constructs the 3D Rotation Matrices R1, R2, R3  (the naming convention is
    // taken from Astrometry), where rotation is performed for the angle "theta"
    // (positive direction is counter-clock-wise as usual) around the X, Y and Z
    // axes, resp. This matrix is actually the INVERSE rotation matrix  (it per-
    // forms the "New-of-Old" co-ords transform), but this is the Astrometrical
    // convention:
    //
    constexpr static Mtx33 MkR1(Angle a_theta)
    {
      double cosTh  = Cos(a_theta);
      double sinTh  = Sin(a_theta);
      Mtx33  R1;
      R1(0,0) = 1.0;    R1(0,1) =  0.0;     R1(0,2) =  0.0;
      R1(1,0) = 0.0;    R1(1,1) =  cosTh;   R1(1,2) =  sinTh;
      R1(2,0) = 0.0;    R1(2,1) = -sinTh;   R1(2,2) =  cosTh;
      return R1;
    }

    constexpr static Mtx33 MkR2(Angle a_theta)
    {
      double cosTh  = Cos(a_theta);
      double sinTh  = Sin(a_theta);
      Mtx33  R2;
      R2(0,0) =  cosTh;   R2(0,1) =  0.0;   R2(0,2) = sinTh;
      R2(1,0) =  0.0;     R2(1,1) =  1.0;   R2(1,2) = 0.0;
      R2(2,0) = -sinTh;   R2(2,1) =  0.0;   R2(2,2) = cosTh;
      return R2;
    }

    constexpr static Mtx33 MkR3(Angle a_theta)
    {
      double cosTh  = Cos(a_theta);
      double sinTh  = Sin(a_theta);
      Mtx33  R3;
      R3(0,0) =  cosTh;   R3(0,1) = sinTh;  R3(0,2) = 0.0;
      R3(1,0) = -sinTh;   R3(1,1) = cosTh;  R3(1,2) = 0.0;
      R3(2,0) =  0.0;     R3(2,1) =  0.0;   R3(2,2) = 1.0;
      return R3;
    }

    //=======================================================================//
    // Low-Level Matrix-Vector Multiplication (Low-Level):                   //
    //=======================================================================//
    // NB:
    // (*) the Matrix is assumed to be of "double"s, but the Vector can be of
    //     any suitable type (typically a "DimQ");
    // (*) we CANNOT use a higher-order API ("Vector3D") here, because the lat-
    //     ter is parameterised by COS, and multiplication by a Rotation Matrix
    //     would change the COS...
    //
    template<typename DQ>
    constexpr void MVMult(DQ const a_vec[3], DQ a_res[3]) const
    {
      assert(a_vec != nullptr && a_res != nullptr);
      for (int i = 0; i < 3; ++i)
      {
        a_res[i] = DQ(0.0);
        for (int j = 0; j < 3; ++j)
          a_res[i] += (*this)(i,j) * a_vec[j];
      }
    }
  };
}
// End namespace SpaceBallistics
