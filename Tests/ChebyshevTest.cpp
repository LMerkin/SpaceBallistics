// vim:ts=2:et
//===========================================================================//
//                             "ChebyshevTest.cpp":                          //
//               Test of Chebyshev Polynomials Implementation                //
//===========================================================================//
#include "SpaceBallistics/Maths/Chebyshev.hpp"
#include <iostream>
#include <cfloat>
#include <ctime>
#include <functional>
#include <random>


int main(int argc, char const* argv[])
{
  using namespace std;
  using namespace SpaceBallistics;
  using namespace DimTypes::Bits::CEMaths;

  //-------------------------------------------------------------------------//
  // Initialisation:                                                         //
  //-------------------------------------------------------------------------//
  // Set the Max Degree if specified in the command line:
  int n = 21;
  if (argc >= 2)
    n = atoi(argv[1]);
  cout << "Degree n = " << n << endl;

  // Initialise the RNG:
  uniform_real_distribution<double> distr(0.0, 1.0);
  ranlux48_base eng(uint64_t(time(nullptr)));
  function<double()> rngr = bind(distr, eng);

  // Testing vars:
  double max_err;
  double res;

  // Generate a random point in [-1.0; +1.0]:
  double x = 2.0 * rngr() - 1.0;
  cout << "Point  x = " << x << endl;
  cout << "Mach.Eps = " << Eps<double> << endl;

  // Vector of random weights in [0.0; +1.0]:
  double  a[n+1];
  for (int j = 0; j <= n; ++j)
    a[j] = rngr();

  //-------------------------------------------------------------------------//
  // TEST 1:                                                                 //
  //-------------------------------------------------------------------------//
  // T_n(zeros_j) == 0.0 indeed:
  //
  cout << "TEST 1: T_n(zeros)    == 0\t\t... ";

  double  zeros[n];
  Chebyshev::Zeros(n, zeros);

  max_err = 0.0;
  for (int j = 0; j < n; ++j)
  {
    res = Abs(Chebyshev::T(n, zeros[j]));

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 2:                                                                 //
  //-------------------------------------------------------------------------//
  // T_0, ..., T_n : one go -vs- one-by-one:
  cout << "TEST 2: T_{0..n}(x)\t\t\t... ";

  double  vals[n+1];
  Chebyshev::Ts(n, x, vals);

  max_err = 0.0;
  for (int j = 0; j <= n; ++j)
  {
    res = Abs(Chebyshev::T(j, x) - vals[j]);

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 3:                                                                 //
  //-------------------------------------------------------------------------//
  // T_n(extrema_j) == +-1.0 indeed:
  //
  cout << "TEST 3: T_n(extrema)  == +-1\t\t... ";

  double  extrema[n+1];
  Chebyshev::Extrema(n, extrema);

  max_err    =  0.0;
  double pm1 =  1.0;

  for (int j = n; j >= 0; --j)
  {
    res = Abs(Chebyshev::T(n, extrema[j]) - pm1);

    if (res > max_err)
      max_err = res;
    pm1 = - pm1;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 4:                                                                 //
  //-------------------------------------------------------------------------//
  // T'_0, ..., T'_n: one go -vs- one-by-one:
  //
  cout << "TEST 4: T'_{0..n}(x)\t\t\t... ";

  double  ders[n+1];
  Chebyshev::DTs(n, x, ders);

  max_err = 0.0;
  for (int j = 0; j <= n; ++j)
  {
    res = Abs(Chebyshev::DT(j, x) - ders[j]);

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 5:                                                                 //
  //-------------------------------------------------------------------------//
  // T'_n(extrema_j) == 0 indeed (internal points only):
  //
  cout << "TEST 5: T'_n(extrema) == 0\t\t... ";
  max_err = 0.0;

  for (int j = 1; j < n; ++j)
  {
    res = Abs(Chebyshev::DT(n, extrema[j]));

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 6:                                                                 //
  //-------------------------------------------------------------------------//
  // T'_{2n+1}(alphas_j) == 0 indeed (j >= 1):
  //
  cout << "TEST 6: T'_{2n+1}(alphas) == 0\t\t... ";

  double  alphas[n+1];
  Chebyshev::Alphas(n, alphas);

  max_err = 0.0;
  int m = 2*n + 1;

  for (int j = 1; j <= n; ++j)
  {
    double xi = 2.0 * alphas[j] - 1.0;
    res = Abs(Chebyshev::DT(m, xi)) / double(m);

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 7:                                                                 //
  //-------------------------------------------------------------------------//
  // T"_0, ..., T"_n: one go -vs- one-by-one:
  //
  cout << "TEST 7: T\"_{0..n}(x)\t\t\t... ";

  double  ddrs[n+1];
  Chebyshev::DDTs(n, x, ddrs);

  max_err = 0.0;
  for (int j = 0; j <= n; ++j)
  {
    res = Abs(Chebyshev::DDT(j, x) - ddrs[j]);

    if (res > max_err)
      max_err = res;
  }
  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 8:                                                                 //
  //-------------------------------------------------------------------------//
  // T"_n(inflects_j) == 0 indeed:
  //
  if (n >= 3)
  {
    cout << "TEST 8: T\"_n(inflects) == 0\t\t... ";

    double  inflects[n-2];
    Chebyshev::Inflects(n, inflects);

    max_err = 0.0;

    for (int j = 0; j < n-2; ++j)
    {
      res = Abs(Chebyshev::DDT(n, inflects[j]));

      if (res > max_err)
      max_err = res;
    }
    cout << "MaxErr = " << max_err << endl;
  }
  else
    cout << "TEST 8: T\"_n(inflects) == 0\t\t... SKIPPED" << endl;

  //-------------------------------------------------------------------------//
  // TEST 9:                                                                 //
  //-------------------------------------------------------------------------//
  cout << "TEST 9: T\"_n(+-1)\t\t\t... ";

  max_err     = Abs(Chebyshev::DDT(n,  1.0) - Chebyshev::DDTpm1(n,  1.0));
  double atm1 = Abs(Chebyshev::DDT(n, -1.0) - Chebyshev::DDTpm1(n, -1.0));
  if (atm1 > max_err)
    max_err = atm1;

  cout << "MaxErr = " << max_err << endl;

  //-------------------------------------------------------------------------//
  // TEST 10:                                                                //
  //-------------------------------------------------------------------------//
  // Sum'_{i=0}^n a_i T_i(x) : 2 ways...
  //
  cout << "TEST 10: Sum'_{i=0}^n a_i T_i(x)\t... ";

  res = - Chebyshev::Sum1T(n, a, x);
  for (int j = 0; j <= n; ++j)
  {
    double c  = (j == 0) ? (0.5 * a[0]) : a[j];
    res += c * Chebyshev::T(j, x);
  }
  cout << "MaxErr = " << Abs(res) << endl;

  //-------------------------------------------------------------------------//
  // TEST 11:                                                                //
  //-------------------------------------------------------------------------//
  // Sum_{i=1}^n a_i T'_i(x): 2 ways...
  //
  cout << "TEST 11: Sum_{i=1}^n a_i T'_i(x)\t... ";

  res = - Chebyshev::SumDT(n, a, x);
  for (int j = 0; j <= n; ++j) // Incl j==0 as well
    res += a[j] * Chebyshev::DT(j, x);

  cout << "MaxErr = " << Abs(res) << endl;

  //-------------------------------------------------------------------------//
  // TEST 12:                                                                //
  //-------------------------------------------------------------------------//
  // Sum_{i=2}^n a_i T"_i(x): 2 ways...
  //
  cout << "TEST 12: Sum_{i=2}^n a_i T\"_i(x)\t... ";

  res = - Chebyshev::SumDDT(n,   a, x);
  for (int j = 0; j <= n; ++j) // Incl j=0, j=1 as well
    res += a[j] * Chebyshev::DDT(j, x);

  cout << "MaxErr = " << Abs(res) << endl;
};
