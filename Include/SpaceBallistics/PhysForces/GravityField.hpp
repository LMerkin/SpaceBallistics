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
#include <format>
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
    assert(r > a_Re);         // Otherwise, divergence may occurs...
    Acc  mainAcc = a_K / r2;  // Main part of the Gravitational Acceleration

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
    constexpr int B = ((N+6)*(N+1))/2;
    assert(gsl_sf_legendre_array_n(size_t(N)) == size_t(B));

    double Ps   [B];
    double DerPs[B];

    // Pre-Compute the Legendre Polynomials and their Derivatives up to order
    // l_max = N:
    int rc =
      gsl_sf_legendre_deriv_array
        (GSL_SF_LEGENDRE_SPHARM, size_t(N), sinPhi, Ps, DerPs);
    assert(rc == 0);

    // Partial Derivatives of the Gravitational Potential wrt x, y, z:
    for (int i = 0; i < 3; ++i)
    {
      double A = double(a_pos[i] / r);
      double B =
        (i == 0)
        ? double(- x * z / r2) :
        (i == 1)
        ? double(- y * z / r2)
        : 1.0 - Sqr(double(z / r));

      double C =
        (i == 0)
        ? double(- r * y / r2xy)  :
        (i == 1)
        ? double(  r * x / r2xy)
        : 0.0;

      // Sum over all Spherical Harmonics:
      double F   = 0.0;
      double irl = Sqr(ir);

      for (int l = 2; l <= N; ++l)
      {
        double mSum = 0.0;
        int    j    = (l*(l+1))/2;
        double Al1  = A * double(l+1);

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
          mSum +=
            (B * P1 - Al1 * P) * (SHC.m_Clm * c + SHC.m_Slm * s) +
            double(m) * C * P  * (SHC.m_Slm * c - SHC.m_Clm * s);
        }
        F   += mSum * irl;
        irl *= ir;
      }
      // NB: In GSL, the Spherical Harmonics are normalised to 1, not to 4*Pi
      // as assumed in the "a_coeffs", so compensate for that:
      F *= SqRt4Pi;

      // And the Main Term:
      a_acc[i] += mainAcc * (F - A);
    }
  }
}
// End namespace SpaceBallistics
