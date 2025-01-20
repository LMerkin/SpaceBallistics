// vim:ts=2:et
//===========================================================================//
//                 "SpaceBallistics/Maths/Dichotomy.hpp":                    //
//                Solving f(x) = 0 by Dichotomy (BiSection)                  //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Types.hpp"
#include <cassert>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Dichotomy" Function:                                                   //
  //=========================================================================//
  template <typename F, typename T>
  constexpr T Dichotomy(F const& a_f, T a_a, T a_b, T a_eps)
  {
    assert(a_a  < a_b);
    auto fa = a_f(a_a);
    auto fb = a_f(a_b);

    constexpr int MaxIters = 1000;
    int i = 0;
    for (;  i < MaxIters; ++i)
    {
      if (UNLIKELY(IsZero(fa)))
        return a_a;
      if (UNLIKELY(IsZero(fb)))
        return a_b;

      // The root must be isolated:
      assert((IsNeg(fa) && IsPos(fb)) || (IsPos(fa) && IsNeg(fb)));

      // Mid-Point:
      T m  = (a_a + a_b) / 2.0;

      // Done?
      if (UNLIKELY(a_b - a_a < a_eps))
        return m;

      // Continue:
      auto fm = a_f(m);
      if (UNLIKELY(IsZero(fm)))
        return m;

      if ((IsNeg(fa) && IsPos(fm)) || (IsPos(fa) && IsNeg(fm)))
      {
        a_b = m;
        fb  = fm;
      }
      else
      {
        assert((IsNeg(fm) && IsPos(fb)) || (IsPos(fm) && IsNeg(fb)));
        a_a = m;
        fa  = fm;
      }
    }
    // If we got here, divergence (XXX???) had occured:
    assert(false);
    return T(NaN<double>);
  }
}
// End namespace SpaceBallistics
