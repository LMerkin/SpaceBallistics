// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/PhysForces/GravityField.h":              //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include "SpaceBallistics/CoOrds/BodyCentricCOSes.h"
#include "SpaceBallistics/PhysForces/BodyData.hpp"
#include <type_traits>
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_sf_legendre.h>

namespace SpaceBallistics
{
  //=========================================================================//
  // "GravityField" Class:                                                   //
  //=========================================================================//
  // Provides the Gravitational Field Model for the given Ellipsoidal Gravita-
  // ting Body (EGB):
  //
  template<Body EGB>
  class GravityField
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // Equatorial Radius of the Body (used in the Gravitational Potential exp-
    // ansion):
    constexpr static LenK Re = BodyData<EGB>::Re;

    // Gravitational Field Constant of the Body:
    constexpr static GMK  K  = BodyData<EGB>::K;

    //=======================================================================//
    // Model Coeffs:                                                         //
    //=======================================================================//
    // The struct for Dimension-Less Spherical Harmonics Coeffs representing the
    // Gravitational Potential, with Geodesy-style normalisation (to 4*Pi):
    //
    struct SpherHarmonicCoeffs
    {
      // Data Flds:
      int    const m_l;      // 2  .. MaxDeg
      int    const m_m;      // 0  .. m_l
      double const m_Clm;    // Coeff at Cos
      double const m_Slm;    // Coeff at Sin
    };

    // The Max Degree and Order of Spherical Harmonics available:
    constexpr static int N = BodyData<EGB>::MaxSpherHarmDegreeAndOrder;

  private:
    // Actual Coeffs:
    static_assert(N == 0 || N >= 2);
    static SpherHarmonicCoeffs const s_coeffs[N == 0 ? 0 : ((N+1)*(N+2))/2];

  public:
    //-----------------------------------------------------------------------//
    // For convenience: Exception thrown on "impact" (actually when r <= Re) //
    //-----------------------------------------------------------------------//
    struct ImpactExn
    {
      Time   const m_t;      // Time
      Len    const m_h;      // Elevation over the Equatorial Radius
      Angle  const m_lambda; // Impact Site Longitude
      Angle  const m_phi;    // Impact Site Latitude
    };


    //=======================================================================//
    // Gravitational Acceleration Computation:                               //
    //=======================================================================//
    // The function ADDS the components of the gravitational  acceleration vec-
    // tor to the output vector "acc",  so the latter must be properly initial-
    // ised (eg zeroed-out) before calling this function.
    // NB: "pos" and "acc" are in the "BodyCRotCOS" (which is embedded in the
    // Body is and rotating with it):
    //
    static void GravAcc
    (
      Time                 a_t,                  // For info only
      PosKVRot<EGB> const& a_pos,
      AccVRot <EGB>*       a_acc,
      int                  a_n          = N,     // Max order used
      bool                 a_zonal_only = false  // Zonal Harmonics only?
    )
    {
      //---------------------------------------------------------------------//
      // Checks:                                                             //
      //---------------------------------------------------------------------//
      static_assert(N >= 0 && IsPos(K) && IsPos(Re));
      assert(a_acc != nullptr);

      if (UNLIKELY(a_n < 0 || a_n == 1))
        throw std::invalid_argument
              ("GravAcc: Invalid Order (must be 0 or >= 2");

      if (UNLIKELY(a_n > N))
        throw std::invalid_argument("GravAcc: Requested Order too high");

      //---------------------------------------------------------------------//
      // The Rectangular and Spherical CoOrds:                               //
      //---------------------------------------------------------------------//
      // (In the "BodyCRotCOS"):
      LenK x       = a_pos[0];
      LenK y       = a_pos[1];
      LenK z       = a_pos[2];
      auto r2xy    = Sqr(x) + Sqr(y);
      auto r2      = r2xy   + Sqr(z);
      LenK r       = SqRt(r2);

      if (UNLIKELY(r <= Re))
      {
        // Inner points are not allowed: Divergence may occur. We treat this as
        // a "surface impact" event,   though it might not be a physical impact
        // yet (we are under the Equatorial Radius, possibly not the local one):
        double const phi     = ASin (double(z/r));
        double const lambda  = IsZero(r2xy) ? 0.0 : ATan2(x, y);

        throw ImpactExn{ a_t, To_Len(r - Re), Angle(lambda), Angle(phi) };
      }

      // If OK: Main part of the Gravitational Acceleration:
      Acc  mainAcc = To_Len_m(K / r2);

      if (a_n == 0)
      {
        // Trivial Case: Spherically-Symmetric Gravitational Field:
        for (int i = 0; i < 3; ++i)
          (*a_acc)[i] = - mainAcc * double(a_pos[i] / r);
        return;
      }
      //---------------------------------------------------------------------//
      // GENERAL CASE: Will sum up the Spherical Harmonics:                  //
      //---------------------------------------------------------------------//
      constexpr double SqRt4Pi = 2.0 * SqRt(Pi<double>);

      // Re / r:
      double const ir = double(Re / r);
      assert(ir < 1.0);

      // The Latitude (via its Sin):
      double const sinPhi  = double(z / r);

      // The Longitude: Undefined if rXY=0, assume Lambda=0 in that case:
      double const lambda  =
        IsZero(r2xy) ? 0.0 : ATan2(x, y);

      // Pre-compute Cos(m*lambda), Sin(m*lambda) for m = 0..a_n:
      assert(2 <= a_n && a_n <= N);

      double cosMLambda[N+1];
      double sinMLambda[N+1];
      for (int m = 0; m <= a_n; ++m)
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
      // l_max = a_n <= N:
      DEBUG_ONLY(int rc =)
        gsl_sf_legendre_deriv_array
          (GSL_SF_LEGENDRE_SPHARM, size_t(a_n), sinPhi, Ps, DerPs) ;
      assert(rc == 0);

      // dr/d{x,y,z}:
      double const A[3]       { double(x/r),   double(y/r),     double(z/r)  };

      // r * d(phi)/d{x,y,z}:
      double const B[3]
        { double(- x*z / r2),   double(- y*z / r2),   1.0 - Sqr(double(z/r)) };

      // r * d(lambda)/d{x,y,z}:
      double const C[3]
        { double(- r*y / r2xy), double(  r*x / r2xy), 0.0 };

      //---------------------------------------------------------------------//
      // Sum over the Spherical Harmonics (l,m):                             //
      //---------------------------------------------------------------------//
      double F[3]  {0.0, 0.0, 0.0};
      double irl = Sqr(ir);

      for (int l = 2; l <= a_n; ++l)
      {
        double mSum[3] {0.0, 0.0, 0.0};
        double l1   =  double(l+1);
        double Al1 [3] {A[0]*l1, A[1]*l1, A[2]*l1};
        int    j    =  (l*(l+1))/2;
        int    maxM =  a_zonal_only ? 0 : l;

        for (int  m = 0; m <= maxM; ++m)
        {
          assert(m <= a_n);
          double c  = cosMLambda[m];
          double s  = sinMLambda[m];

          // Index for accessing Legendre Polys and Model Coeffs:
          int    jm = j + m;
          assert(jm < ((a_n+1)*(a_n+2))/2);

          SpherHarmonicCoeffs SHC = s_coeffs[jm];
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
      //---------------------------------------------------------------------//
      // Finally:                                                            //
      //---------------------------------------------------------------------//
      for (int i = 0; i < 3; ++i)
      {
        // NB: In GSL, the Spherical Harmonics are normalised to 1, not to 4*Pi
        // as assumed in the "s_coeffs", so compensate for that:
        F[i] *= SqRt4Pi;

        // And the Main Term:
        (*a_acc)[i] += mainAcc * (F[i] - A[i]);
      }
    }
  };
}
// End namespace SpaceBallistics
