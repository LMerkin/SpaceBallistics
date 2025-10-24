// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Maths/RKF5.hpp":                      //
//    DimTypes-Enabled Runge-Kutta-Fehlberg ODE Integrator (of Order 4/5)    //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Maths/LinAlgT.hpp"
#include <stdexcept>
#include <ostream>

namespace SpaceBallistics
{
  //=========================================================================//
  // "RKF5":                                                                 //
  //=========================================================================//
  // ODE solution using the Runge-Kutta-Fehlberg method,  with automatic step
  // control.
  // Returns the last "time" instant which the output "a_x" corresponds to; if
  // it equals "a_tf", integration is successfully completed:
  //
  template<typename X, typename T, typename RHS, typename CB>
  T RKF5
  (
    X*            a_x,       // Initial condition, replaced by the solution
    T             a_t0,      // Initial independent variable (eg "time")
    T             a_tf,      // Final "time"
    RHS const&    a_rhs,     // RHS function:   (x,t) -> dx/dt
    T             a_tau0,    // Initial "time" step
    T             a_max_tau, // Max     "time" step
    double        a_eps,     // Max allowed absolute error
    CB*           a_cb,      // User Call-Back: (&x,t,tau)->ContFlag (NULL OK)
    std::ostream* a_log      // Stream for logging (may be NULL)
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    if (a_x == nullptr ||
       (a_tf > a_t0 && !IsPos(a_tau0)) || (a_tf < a_t0 &&  IsPos(a_tau0)) ||
       Abs(a_tau0)  >  Abs(a_max_tau))
      throw std::invalid_argument("RKF5: Invalid param(s)");

    if (a_tf == a_t0)
      return a_t0;     // Nothing to do

    DEBUG_ONLY(bool isFwd = a_tf > a_t0;)

    //-----------------------------------------------------------------------//
    // Outer Loop: Time Marshaling:                                          //
    //-----------------------------------------------------------------------//
    T t   = a_t0;
    T tau = a_tau0;

    while (true)
    try    // Guard against any exceptions, eg in the RHS evaluation
    {
      // Invoke the Call-Back if present, stop immediately if it returns
      // "false". If there is no Call-Back, just continue.
      // NB: The Call-Back MAY (in some cases) modify "a_x":
      //
      if (a_cb != nullptr && !(*a_cb)(a_x, t, tau))
        return t;     // We may or may not have reached the end

      if (t == a_tf)
        return t;     // All Done!

      // The next "t" to be reached in one RKF5 step:
      T tn = t + tau;

      // But if we can now reach "a_tf" in one step, do so:
      if (UNLIKELY(Abs(a_tf - t) <= 1.01 * Abs(tau)))
      {
        tn  = a_tf;
        tau = a_tf - t;
      }
      assert((isFwd && tn <= a_tf) || (!isFwd && tn >= a_tf));

      //---------------------------------------------------------------------//
      // Inner Loop: t -> tn, with possible step reductions:                 //
      //---------------------------------------------------------------------//
      // XXX: Allow no more than 10 step reductions; if after that we still
      // cannot meet the accuracy requirements, something is deeply wrong!
      // NB: the step can only be reduced, never enlarged:
      //
      auto f0 = a_rhs(*a_x, t, tau);
      bool ok = false;

      for (int nsr = 0; nsr < 10; ++nsr)
      {
        X    x1   = Add (*a_x, 0.25 * tau, f0);
        T    t1   = t + 0.25 * tau;
        auto f1   = a_rhs(x1, t1, tau);

        auto b2   = Mult(     0.09375, f0);
        b2        = Add (b2,  0.28125, f1);
        X    x2   = Add (*a_x,  tau,   b2);
        T    t2   = t + 0.375 * tau;
        auto f2   = a_rhs(x2, t2, tau);

        auto b3   = Mult(      1932.0 / 2197.0, f0);
        b3        = Add (b3,  -7200.0 / 2197.0, f1);
        b3        = Add (b3,  7296.0  / 2197.0, f2);
        X    x3   = Add (*a_x, tau, b3);
        T    t3   = t + (12.0 / 13.0) * tau;
        auto f3   = a_rhs(x3, t3, tau);

        auto b4   = Mult(      439.0  /  216.0, f0);
        b4        = Add (b4,  -8.0,             f1);
        b4        = Add (b4,  3680.0  /  513.0, f2);
        b4        = Add (b4,  -845.0  / 4104.0, f3);
        X    x4   = Add (*a_x, tau, b4);
        // NB: t4=tn:
        auto f4   = a_rhs(x4, tn, tau);

        auto b5   = Mult(     -8.0    /   27.0, f0);
        b5        = Add (b5,   2.0,             f1);
        b5        = Add (b5,  -3544.0 / 2565.0, f2);
        b5        = Add (b5,   1859.0 / 4104.0, f3);
        b5        = Add (b5,  -0.275,           f4);
        X    x5   = Add (*a_x, tau, b5);
        T    t5   = t + 0.5  * tau;
        auto f5   = a_rhs(x5, t5, tau);

        // Error estimate (NB: "f1" is not used here, this is correct):
        auto ef   = Mult(        1.0 / 360.0,   f0);
        ef        = Add (ef,  -128.0 / 4275.0,  f2);
        ef        = Add (ef, -2197.0 / 75240.0, f3);
        ef        = Add (ef,    0.02,           f4);
        ef        = Add (ef,     2.0 / 55.0,    f5);
        X    ex   = Mult(tau, ef);

        // NB: Max Abs Error (with "tau" factor), as a plain "double", since
        // the actual element type on which this max occurs, cannot be known
        // in advance:
        double maxAbsErr  = MaxAbsRelError(ex, x5);
        assert(maxAbsErr >= 0.0);

        // XXX: If eps <= 0, do not use error / step control at all:
        if (a_eps <= 0.0 || maxAbsErr < a_eps)
        {
          // OK, the error bounds are satisfied, update "a_x". NB: "f1" is not
          // used, this is correct:
          auto b6 = Mult(       16.0 /   135.0, f0);
          b6      = Add (b6,  6656.0 / 12825.0, f2);
          b6      = Add (b6, 28561.0 / 56430.0, f3);
          b6      = Add (b6,   -0.18,           f4);
          b6      = Add (b6,     2.0 / 55.0,    f5);

          // This step is done:
          *a_x    = Add(*a_x, tau, b6);
          t       = tn;
          ok      = true;

          // Furthermore, we can enlarge "tau", but not beyond "a_max_tau":
          if (a_eps > 0.0 && tau <= 0.5 * a_max_tau)
          {
            T tau1 =
              Min(Abs(tau) * std::pow(a_eps / maxAbsErr, 0.2), Abs(a_max_tau));
            assert(IsPos(tau1));

            // Increase "tau" by 2-factors, up to "tau1":
            while (true)
            {
              T tau2 = 2.0  * tau;
              if (Abs(tau2) > tau1)
                break;
              tau = tau2;
            }
          }
          // The inner loop (one RKF step) is done!
          break;
        }
        else
        {
          // Error bound "eps" exceeded; we need to reduce the step and try
          // again. New step estimate:
          assert(a_eps > 0.0);
          T tau1 = 0.9 * Abs(tau) * std::pow(a_eps / maxAbsErr, 0.2);
          assert(IsPos(tau1));

          // However, we will not use "tau1" directly; better reduce the origi-
          // nal "tau" by powers of 2 until it becomes <= "tau1":
          while (Abs(tau) > tau1)
            tau *= 0.5;

          // And go for the next inner iteration with the new "tau"...
        }
      }
      // End of the inner (single RKF5 step) loop.
      // If all step reductions have been done, but to no avail, throw an
      // exception:
      if (!ok)
        throw std::logic_error("RKF5: Too many step reductions");
    }
    catch (std::exception const& exn)
    {
      if (a_log != nullptr)
        *a_log  << "# RKF5: Exception: " << exn.what() << std::endl;
      throw;
    }
    // End of (exception-protected) "t" loop
    // This point is unreachable:
    __builtin_unreachable();
  }
}
// End namespace SpaceBallistics
