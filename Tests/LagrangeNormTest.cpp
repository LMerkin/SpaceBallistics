// vim:ts=2:et
//===========================================================================//
//                       "Tests/LagrangeNormTest.cpp":                       //
//  Test the Normalisation Coeffs of Associated Lagrange Polynomials in GSL  //
//===========================================================================//
#include "DimTypes/Bits/Macros.h"
#include <gsl/gsl_sf_legendre.h>
#include <iostream>
#include <cassert>
#include <cmath>

int main()
{
  using namespace std;

  // Allocate the result space:
  constexpr int N  = 10;
  constexpr int NP = ((N+6)*(N+1))/2;
  assert(gsl_sf_legendre_array_n(size_t(N)) == size_t(NP));
  double Ps    [NP];

  // Compute the Associated Legendre Polynomials @ x=1 for dgrees 0..N:
  DEBUG_ONLY(int rc =)
    gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM, size_t(N), 1.0, Ps);
  assert(rc == 0);

  // The normalisation coeffs for P(l,0) (n=0..N) should be
  // SqRt((2*l+1)/(4*Pi)): Check this:
  //
  for (int l = 0; l <= N; ++l)
  {
    double val = Ps[(l*(l+1))/2];
    cout << "l=" << l << ": diff=" << (val - sqrt(double(2*l+1)/(4.0*M_PI)))
         << endl;
  }
  return 0;
}
