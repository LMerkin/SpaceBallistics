// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Maths/Chebyshev.hpp":                 //
//      Computation of Chebyshev polynomials, derivatives and series         //
//===========================================================================//
#pragma  once

#include <DimTypes/Bits/CEMaths.hpp>
#include <DimTypes/Bits/Macros.h>
#include <type_traits>
#include <cassert>

namespace SpaceBallistics::Chebyshev
{
  using namespace DimTypes::Bits::CEMaths;

  //=========================================================================//
  //  T(n, x):                                                               //
  //=========================================================================//
  //  Value of the Chebyshev polynomial of order "n" (n >= 0) at "x":
  //
  template<typename F>
  constexpr F T(int n, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && Abs(x) <= F(1.0));

    switch (n)
    {
      case  0: return F(1.0);
      case  1: return x;
      default:
      {
        F T0    = F(1.0);
        F T1    = x;
        F two_x = F(2.0) * x;

        for (int i = 2; ; ++i)
        {
          F T2 = two_x * T1 - T0;

          if (i == n)
            return T2;

          T0 = T1;
          T1 = T2;
        }
      }
    }
  }

  //=========================================================================//
  //  Ts(n, x):                                                              //
  //=========================================================================//
  //  T_0(x), ..., T_n(x) evaluated at the same point "x"; results
  //  stored in "res"; len(res) == n+1 (must have sufficient size):
  //
  template<typename F>
  constexpr void Ts(int n, F x, F* res)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && res != nullptr && Abs(x) <= F(1.0));

    F two_x = F(2.0) * x;

    for (int i = 0; i <= int(n); ++i)
      switch (i)
      {
        case 0:
          res[i] = F(1.0);
          break;
        case 1:
          res[i] = x;
          break;
        default:
          res[i] = two_x * res[i-1] - res[i-2];
      }
  }

  //=========================================================================//
  //  Sum1T(n, a, x):                                                        //
  //=========================================================================//
  //  Efficient evaluation of Sum'_{i=0}^n a_i T_i(x), len(a) == n+1;
  //  NB: "L" may be any linear space over "F":
  //
  template<typename F, typename L>
  constexpr L Sum1T(int n, L const* a, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && a != nullptr && Abs(x) <= F(1.0));

    L B1    = L(0.0);
    L B2    = L(0.0);
    F two_x = F(2.0) * x;

    for (int i = n; ; --i)
    {
      L B0 = two_x * B1 - B2 + a[i];

      if (i == 0)
        return F(0.5) * (B0 - B2);

      B2 = B1;
      B1 = B0;
    }
  }

  //=========================================================================//
  //  DT(n, x):                                                              //
  //=========================================================================//
  //  dT_n/dx at x, n >= 0:
  //
  template<typename  F>
  constexpr F DT(int n, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && Abs(x) <= F(1.0));

    switch (n)
    {
      // T_0(x) = 1, T_1(x) = x, T_2(x) = 2 x^2 - 1, so:
      case 0:  return  F(0.0);
      case 1:  return  F(1.0);
      case 2:  return  F(4.0) * x;
      default:
      // Generic case: n >= 3, d T_n/dx = n U_{n-1}, where
      // U(x) is the Chebyshev polynomial of 2nd kind:
      {
        F U0    = F(1.0);
        F two_x = F(2.0) * x;
        F U1    = two_x;

        for (int i = 3; ; ++i)
        {
          F U2 = two_x * U1 - U0;

          if (i == n)
            return F(n) * U2;

          U0 = U1;
          U1 = U2;
        }
      }
    }
  }

  //=========================================================================//
  // DTs(n, x):                                                              //
  //=========================================================================//
  // dT_i/dx at x, for all i=0..n; res[0] is always set to 0;
  // pre-cond: len(res) == n+1:
  //
  template<typename F>
  constexpr void DTs(int n, F x, F* res)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && res != nullptr && Abs(x) <= F(1.0));

    F U0    = F(1.0);
    F two_x = F(2.0) * x;
    F U1    = two_x;

    for (int i = 0; i <= n; ++i)
    {
      // T_0(x) = 1, T_1(x) = x, T_2(x) = 2 x^2 - 1, so:
      switch (i)
      {
        case 0:
          res[i] = F(0.0);
          break;
        case 1:
          res[i] = F(1.0);
          break;
        case 2:
          res[i] = F(4.0) * x;
          break;
        default:
        {
          // Generic case: i >= 3, D t_i/dx = i U_{i-1}, where
          // U(x) is the Chebyshev polynomial of 2nd kind:
          F  U2  = two_x * U1 - U0;
          res[i] = F(i)  * U2;
          U0     = U1;
          U1     = U2;
        }
      }
    }
  }

  //=========================================================================//
  //  SumDT(n, a, x):                                                        //
  //=========================================================================//
  //  Efficient evaluation of Sum_{i=1}^{n} a_i dT_i/dx, len(a) == n+1,
  //  a[0] is irrelevant, so in this case Sum' = Sum;
  //  Again, "L" may be any linear space over "F":
  //
  template<typename F,  typename L>
  constexpr L SumDT(int n, L const* a, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && a != nullptr && Abs(x) <= F(1.0));

    // We need to compute: Sum_{i=1}^n i a_i U_{i-1}(x), where U_i(x) are
    // Chebyshev polynomials of 2nd kind; a_0 is unused.
    // The recursive procedure is similar to that of "Sum1T",    only the
    // final return value is different:
    //
    if (n == 0)
      return L(0.0);

    L B1    = L(0.0);
    L B2    = L(0.0);
    F two_x = F(2.0) * x;

    for (int i = n; ; --i)
    {
      L B0 = two_x * B1 - B2 + F(i) * a[i];

      if (i == 1)
        return B0;

      B2 = B1;
      B1 = B0;
    }
  }

  //=========================================================================//
  //  DDT(n, x):                                                             //
  //=========================================================================//
  //  d^2 T_n / dx^2 at "x", n >= 0:
  //
  template<typename F>
  constexpr F DDT(int n, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && Abs(x) <= F(1.0));

    if (n <= 1)
      return F(0.0);

    // Generic case. The formula is:
    //
    // Sum'_{i >= 0, i <= n-2, n-i even} (n-i) n (n+i) T_i(x)
    //
    // where "'" means: coeff of T_0 is multiplied by F(0.5) :
    //
    F B1    = F(0.0);
    F B2    = F(0.0);
    F two_x = F(2.0) * x;
    bool   even  = true;

    for (int i = n-2; ; --i)
    {
      F B0 = two_x * B1 - B2;
      if (even)
        B0 += F((n-i) * n * (n+i));

      if (i == 0)
        return F(0.5) * (B0 - B2);

      B2 = B1;
      B1 = B0;
      even = !even;
    }
  }

  //=========================================================================//
  //  DDTpm1(n, x):                                                          //
  //=========================================================================//
  //  Special case of the above: "x" must be +-1 :
  //
  template<typename F>
  constexpr F DDTpm1(int n, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && Abs(x) == F(1.0));

    if (n <= 1)
      return F(0.0);

    // We sum up the series of either even- or odd-indiced T_j (depending on
    // "n"). Their value at +1 is always +1, at -1 -- depends on "n":
    //
    F sign = (x == F(1.0) || n % 2 == 0) ? F(1.0) : -F(1.0);
    F res  = F(0.0);

    for (int r = n-2; r >= 0; r -= 2)
      if (r == 0)
        res += F(0.5)  * F(n * n * n);
      else
        res += sign * F((n-r) * n * (n+r));
    return res;
  }

  //=========================================================================//
  //  DDTs(n, x):                                                            //
  //=========================================================================//
  // d^2 T_i / dx^2 at x,  for all  i=0..n;
  // res[0] and res[1] are always set to 0;
  // pre-cond: len(res) == n+1:
  //
  template<typename F>
  constexpr void DDTs(int n, F x, F* res)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && Abs(x) <= F(1.0));

    // For i >= 2, the formula is:
    //
    // d^2 T_n / dx^2 = Sum'_{i >= 0, i <= n-2, n-i even} (n-i) n (n+i) T_i(x)
    //
    // where "'" means: coeff of T_0 is multiplied by F(0.5) :
    //
    // Thus: store the Chebyshev polynomials in "res". We only need degrees up
    // to (n-2):
    Ts(n-2, x, res);

    // XXX: still could not avoid a F loop...
    for (; n >= 0; --n)
    {
      // The following loop will over-write "res" from top down:
      res[n] = F(0.0);
      for (int i = n-2; i >= 0; i -= 2)
      {
        F c = F((n-i)*(n+i));
        if (i == 0)
          c /= F(2.0);
        res[n] += c * res[i];
      }
      res[n] *= F(n);
    }
  }

  //=========================================================================//
  //  SumDDT(a, x):                                                          //
  //=========================================================================//
  //  Efficient evaluation of Sum'_{i=2}^{n} a_i d^2 T_i/dx, len(a) == n+1,
  //  a[0] and a[1] are actually ignored, so in this case Sum' = Sum;
  //  Once again,  "A" may be any linear space over "F":
  //
  template<typename F, typename  L>
  constexpr static  L SumDDT(int n, L const* a, F x)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && a != nullptr && Abs(x) <= F(1.0));

    if (n <= 1)
      return F(0.0);

    // Generic case. The formula is:  Sum_{r=1}^{n-1} c_r T'_r(x),
    // where
    // c_r = k_r Sum_{i=r+1, i-r odd }^n a_i i ,
    //       k_0 = 1, k_r = 2 for r >= 1 .
    // We can therefore re-use the "SumDT" algorithm:
    //
    L B1       = L(0.0);
    L B2       = L(0.0);
    L cO       = L(0.0);
    L cE       = L(0.0);
    F    two_x = F(2.0) * x;
    bool odd   = true;

    for (int i = n-1; ; --i)
    {
      L& c  = odd ? cO : cE;
      c   += F(2*(i+1)) * a[i+1];

      L B0 = two_x * B1 - B2 + F(i) * c;

      if (i == 1)
        return B0;

      B2  = B1;
      B1  = B0;
      odd = !odd;
    }
  }

  //=========================================================================//
  //  Zeros(n):                                                              //
  //=========================================================================//
  //  Fills in the "res" vector with "n" real zeros of T_n in [-1.0; +1.0].
  //  Pre-condition: len(res) == n:
  //
  template<typename F>
  constexpr void Zeros(int n, F* res)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && res != nullptr);

    if (n == 0)
      return; // No zeros, but not an error!

    F pi_n = Pi<F> / F(n);

    for (int i = 0; i < n; ++i)
      res[i] = Cos((F(n-i) - F(0.5)) * pi_n);
  }

  //=========================================================================//
  //  Alphas(n):                                                             //
  //=========================================================================//
  //  Fills in "res" with (n+1) odd zeros of the Translated 2nd-kind Chebyshev
  //  polynomial Us_{2n}(x) of degree (2n);  required for the Everhart-Sorokin
  //  integrator. Pre-condition: len(res) == n+1:
  //
  template<typename F>
  constexpr void Alphas(int n,  F* res)
  {
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 0 && res != nullptr);

    F pi_2n1 = Pi<F> / F(2*n+1);
    int    m = 1;

    res[0] = F(0.0);
    for (int j = 1; j <= n; ++j)
    {
      res[j] = F(0.5) * (F(1.0) + Cos(F(m) * pi_2n1));

      assert(res[j] > F(0.0) && res[j] < F(1.0));
      m += 2;
    }
  }

  //=========================================================================//
  //  Extrema(n):                                                            //
  //=========================================================================//
  //  Fills in the "res" vector with (n+1) extremum points of T_n in
  //  [-F(1.0); +F(1.0)], n >= 1.  Raises exception for n <= 0. The end points
  //  (-1, +1) are included in the result: although T'(x) does not vanish
  //  there, T(x) still takes  its min/max  value  (+-1) there -- same as
  //  in internal extremal points. Pre-condition:  len(res) == n+1:
  //
  template<typename F>
  constexpr void Extrema(int n, F* res)
  {
    // There are (n-1) internal extremal points and 2 end points, so (n+1).
    // We disallowed n==0, as in that case T_n(x)==1 identically, and it's
    // unclear which value to assign to the single point (n+1 == 1).    In
    // this case, return an empty vector:
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 1 && res != nullptr);

    F pi_n = Pi<F> / F(n);

    for (int i = 0; i <= n ; ++i)
      res[i] = Cos(F(n-i) * pi_n);
  }

  //=========================================================================//
  //  Inflects(n):                                                           //
  //=========================================================================//
  //  Fills in the "res" vector   with (n-2) inflection points of T_n  in
  //  [-F(1.0); +F(1.0)], i.e. zeros of T"(x). Pre-condition: len(res) == n-2:
  //
  template<typename F>
  constexpr void Inflects(int n, F* res)
  {
    // The problem is more complex here than for Zeros and Extrema, since
    // Inflection Points can only be computed numerically:
    static_assert(std::is_floating_point_v<F>);
    assert(n >= 3 && res != nullptr);

    // "n" even: n = 2*k,   k >= 2;
    // "n" odd : n = 2*k+1, k >= 1;  in any case:
    int k = n / 2;
    assert(k >= 1);

    //  Construct the inflection points {x = cos z}, z in [0..Pi].
    //  There are (k-1) "z" points in (0..Pi/2) and (Pi/2..Pi);
    //  for "n" odd, there is also a point at Pi/2 itself (x==0):
    //
    F dn        = F(n);
    F const tol = (n <= 70) ? F(1e+5 * Eps<F>)
                            : F(1e+8 * Eps<F>);
    // Need such tolerance for large "n"

    for (int m = 1; m < k ; ++m)   // m = [1..k-1]
    {
      // Find the point localised in Pi/(2*n)* (2*m..2*m+1).
      // It is a root of equation n*tan(z) = tan(n*z) . The
      // initial approximation is obtained by the following
      // iterative step:
      F z = F(0.5)   * Pi<F>   * F(2*m+1)      / dn;
      z   = (ATan(dn * Tan(z)) + F(m) * Pi<F>) / dn;

      // Now the iterations:
      // z_k  = z_{k-1} - g(z_k) / g'(z_k), where
      // g(z) = n*tan(z) - tan(n*z) ;
      // g(z) > 0 during the iterative process:
      //
      DEBUG_ONLY(bool ok = false;)

      for(int i = 0; i < 100; ++i)
      {
        F nz     = dn * z;
        F tan_nz = Tan(nz);
        F tan_z  = Tan(z);
        F g      = tan_nz - dn * tan_z;

        // Theoretically, g -> +0 in this process, but due to rounding errors,
        // it could also become negative:
        if (Abs(g) < tol)
        {
          DEBUG_ONLY(ok = true;)
          break;      // "z0" is a root required
        }
        else
        if (g < F(0.0))
          break;      // Jumped into g < 0, would now diverge!

        // Otherwise: OK: compute the next point:
        F dg     = dn * (tan_nz - tan_z) * (tan_nz + tan_z);
        assert(dg > F(0.0));

        g /= dg;
        // "g" is the decrement of "z"; if it is smaller than the same "tol",
        // we exit as well:
        if (g < tol)
        {
          DEBUG_ONLY(ok = true;)
          break;
        }
        z -= g;
      }
      // Iterative solver finished, hopefully successfully:
      assert(ok);

      // Root successfully computed -- place it, and its symmetric
      // counter-part, in the resulting vector, as "x". Increasing
      // values of "m" correspond to increasing NEGATIVE values of
      // "x":
      F x = Cos(z);

      int l = m-1;
      int r = n-2-m;
      assert(l < r);

      res[l] = - x;
      res[r] =   x;
    }
    // Mid-point for odd "n", corresp to m == k, l == r == k-1 :
    if (n % 2)
      res[k-1] = F(0.0);
  }
}
// End namespace SpaceBallistics::Chebyshev
