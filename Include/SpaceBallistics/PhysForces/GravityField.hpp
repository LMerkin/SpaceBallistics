// vim:ts=2:et
//===========================================================================//
//                  "SpaceBallistics/PhysForces/GravityField.h":             //
//===========================================================================//
#pragma once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/GeoCentricRotatingCOS.h"
#include "SpaceBallistics/CoOrds/SelenoCentricRotatingCOS.h"
#include <type_traits>
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_sf_legendre.h>

namespace SpaceBallistics
{
  //-------------------------------------------------------------------------//
  // Model Data:                                                             //
  //-------------------------------------------------------------------------//
  // The struct for Dimension-Less Spherical Harmonics Coeffs  representing  a
  // Gravitational Potential. Geodesy-style normalisation (to 4*Pi) is assumed.
  // Obviously, they depend on the Gravitating Body specified via its Rotating
  // Co-Ords System (eg
  //
  template<typename RCOS>
  struct SpherHarmonicCoeffs
  {
    // XXX: Currently, we have GeoCentric or SelenoCentric Rotating COSes:
    static_assert(std::is_same_v<RCOS, GeoCentricRotatingCOS>   ||
                  std::is_same_v<RCOS, SelenoCentricRotatingCOS>);

    // Data Flds:
    int    const m_l;     // 2  .. MaxDeg
    int    const m_m;     // 0  .. m_l
    double const m_Clm;   // Coeff at Cos
    double const m_Slm;   // Coeff at Sin
  };

  //-------------------------------------------------------------------------//
  // Gravitational Acceleration Computation:                                 //
  //-------------------------------------------------------------------------//
  // The function ADDS the components of the gravitational acceleration vector
  // (which is due to the gravitation field of the Moon)  to the output vector
  // "acc", so the latter must be properly initialised (eg zeroed-out)  before
  // calling this function:
  //
  template<int N, typename RCOS>
  void GravAcc
  (
    GM                              a_K,       // Grav Fld Const
    Len                             a_Re,      // Equator. Radius
    SpherHarmonicCoeffs<RCOS> const (&a_coeffs)[((N+1)*(N+2))/2],
    PosV<RCOS> const&               a_pos,
    AccV<RCOS>*                     a_acc
  )
  {
    static_assert(N == 0 || N >= 2);           // No models of order 1...
    assert(IsPos(a_K) && IsPos(a_Re) && a_acc != nullptr);

    constexpr double SqRt4Pi = 2.0 * SqRt(Pi<double>);

    // The the CoOrds:
    Len  x       = a_pos[0];
    Len  y       = a_pos[1];
    Len  z       = a_pos[2];
    Len2 r2xy    = Sqr(x) + Sqr(y);
    Len2 r2      = r2xy   + Sqr(z);
    Len  r       = SqRt(r2);

    if (UNLIKELY(r <= a_Re))
      // Inner points are not allowed: Divergence may occur:
      throw std::runtime_error("GravAcc: Under the surface");

    // Main part of the Gravitational Acceleration:
    Acc  mainAcc = a_K / r2;

    // Re / r:
    double ir = double(a_Re / r);
    assert(ir < 1.0);

    // The Latitude (via its Sin):
    double sinPhi  = double(z / r);

    // The Longitude: Undefined if rXY=0, assume Lambda=0 in that case:
    double lambda  =
      IsZero(r2xy) ? 0.0 : std::atan2(x.Magnitude(), y.Magnitude());

    // Pre-compute Cos(m*lambda), Sin(m*lambda) for m = 0..N:
    double cosMLambda[N+1];
    double sinMLambda[N+1];
    for (int m = 0; m <= N; ++m)
    {
      double ml = double(m) * lambda;
      cosMLambda[m] = Cos(ml);
      sinMLambda[m] = Sin(ml);
    }
    // Buffers for Legendre Polynomials and Derivatives computation, as per
    // GSL requirements. XXX: Is stack overflow possible here?
    constexpr int NP = ((N+6)*(N+1))/2;
    assert(gsl_sf_legendre_array_n(size_t(N)) == size_t(NP));
    double Ps    [NP];
    double DerPs [NP];

    // Pre-Compute the Legendre Polynomials and their Derivatives up to order
    // l_max = N:
    DEBUG_ONLY(int rc =)
      gsl_sf_legendre_deriv_array
        (GSL_SF_LEGENDRE_SPHARM, size_t(N), sinPhi, Ps, DerPs) ;
    assert(rc == 0);

    // dr/d{x,y,z}:
    double const A[3]       { double(x/r),   double(y/r),     double(z/r)  };

    // r * d(phi)/d{x,y,z}:
    double const B[3]
      { double(- x*z / r2),   double(- y*z / r2),   1.0 - Sqr(double(z/r)) };

    // r * d(lambda)/d{x,y,z}:
    double const C[3]
      { double(- r*y / r2xy), double(  r*x / r2xy), 0.0 };

    // Sum over all Spherical Harmonics:
    double F[3]  {0.0, 0.0, 0.0};
    double irl = Sqr(ir);

    for (int l = 2; l <= N; ++l)
    {
      double mSum[3] {0.0, 0.0, 0.0};
      double l1   =  double(l+1);
      double Al1 [3] {A[0]*l1, A[1]*l1, A[2]*l1};
      int    j    =  (l*(l+1))/2;

      for (int  m = 0; m <= l; ++m)
      {
        assert(m <= N);
        double c  = cosMLambda[m];
        double s  = sinMLambda[m];

        // Index for accessing Legendre Polys and Model Coeffs:
        int    jm = j + m;
        assert(jm < ((N+1)*(N+2))/2);

        SpherHarmonicCoeffs SHC = a_coeffs[jm];
        assert(SHC.m_l == l && SHC.m_m == m);
        double P  = Ps   [jm];
        double P1 = DerPs[jm];

        for (int i = 0; i < 3; ++i)
          mSum[i] +=
            (B[i] * P1 - Al1[i] * P) * (SHC.m_Clm * c + SHC.m_Slm * s) +
            double(m)  * C[i]   * P  * (SHC.m_Slm * c - SHC.m_Clm * s);
      }
      for (int i = 0; i < 3; ++i)
        F[i] += mSum[i] * irl;
      irl *= ir;
    }
    // Finally:
    for (int i = 0; i < 3; ++i)
    {
      // NB: In GSL, the Spherical Harmonics are normalised to 1, not to 4*Pi
      // as assumed in the "a_coeffs", so compensate for that:
      F[i] *= SqRt4Pi;

      // And the Main Term:
      (*a_acc)[i] += mainAcc * (F[i] - A[i]);
    }
  }
}
// End namespace SpaceBallistics
