// vim:ts=2:et
//===========================================================================//
//                    "SpaceBallistics/Maths/RKF5.hpp":                      //
//    DimTypes-Enabled Runge-Kutta-Fehlberg ODE Integrator (of Order 4/5)    //
//===========================================================================//
#pragma  once
#include "SpaceBallistics/Maths/LinAlgT.hpp"
#include <stdexcept>

namespace SpaceBallistics
{
  //=========================================================================//
  // "RKF5":                                                                 //
  //=========================================================================//
  // ODE solution using the Runge-Kutta-Fehlberg method, with automatic step
  // control:
  //
  template<typename X, typename T, typename RHS, typename CB>
  void RKF5
  (
    X*         a_x,    // Initial condition, replaced by the solution
    T          a_t0,   // Initial independent variable (eg "time")
    T          a_tf,   // Final "time"
    RHS const& a_rhs,  // RHS function:   (x,t) -> dx/dt
    T          a_tau,  // Initial "time" step
    double     a_eps,  // Max allowed absolute error
    CB  const& a_cb    // User Call-Back: (x,t) -> continue_flag
  )
  {
    //-----------------------------------------------------------------------//
    // Checks:                                                               //
    //-----------------------------------------------------------------------//
    if (a_x == nullptr ||
       (a_tf > a_t0 && !IsPos(a_tau)) || (a_tf < a_t0 &&  IsPos(a_tau)))
      throw std::invalid_argument("RKF5: Invalid param(s)");

    if (a_tf == a_t0)
      return;   // Nothing to do

    bool isFwd = a_tf > a_t0;

    //-----------------------------------------------------------------------//
    // Outer Loop: Time Marshaling:                                          //
    //-----------------------------------------------------------------------//
    T      t      = a_t0;
    T      tau    = a_tau;
    double tauMag = double(tau / T(1.0));

    while (true)
    {
      // Invoke the Call-Back. Stop immediately if it returns "false":
      if (!a_cb(*a_x, t))
        break;

      if (t == a_tf)
        break;   // All Done!

      // The next "t" to be reached in one RKF5 step:
      T tn =
        std::fabs(tau) >= 0.99 * std::fabs(a_tf - t)
        ? a_tf      // We can reach "a_tf" in one step!
        : t + tau;  // Just an intermediate step
      assert((isFwd && tn <= a_tf) || (!isFwd && tn >= a_tf));

      //---------------------------------------------------------------------//
      // Inner Loop: t -> tn, with possible step reductions:                 //
      //---------------------------------------------------------------------//
      // XXX: Allow no more than 10 step reductions; if after that we still
      // cannot meet the accuracy requirements, something is deeply wrong!
      // NB: the step can only be reduced, never enlarged:
      //
      auto f0 = a_rhs(*a_x, t);
      bool ok = false;

      for (int nsr = 0; nsr < 10; ++nsr)
      {
        X    x1   = Add (*a_x, 0.25 * tau, f0);
        T    t1   = t + 0.25 * tau;
        auto f1   = a_rhs(x1, t1);

        auto b2   = Mult(     0.09375, f0);
        b2        = Add (b2,  0.28125, f1);
        X    x2   = Add (*a_x,  tau,   b2);
        T    t2   = t + 0.375 * tau;
        auto f2   = a_rhs(x2, t2);

        auto b3   = Mult(      1932.0 / 2197.0, f0);
        b3        = Add (b3,  -7200.0 / 2197.0, f1);
        b3        = Add (b3,  7296.0  / 2197.0, f2);
        X    x3   = Add (*a_x, tau, b3);
        T    t3   = t + (12.0 / 13.0) * tau;
        auto f3   = a_rhs(x3,  t3);

        auto b4   = Mult(      439.0  /  216.0, f0);
        b4        = Add (b4,  -8.0,             f1);
        b4        = Add (b4,  3680.0  /  513.0, f2);
        b4        = Add (b4,  -845.0  / 4104.0, f3);
        X    x4   = Add (*a_x, tau, b4);
        // NB: t4=tn:
        auto f4   = a_rhs(x4,  tn);

        auto b5   = Mult(     -8.0    /   27.0, f0);
        b5        = Add (b5,   2.0,             f1);
        b5        = Add (b5,  -3544.0 / 2565.0, f2);
        b5        = Add (b5,   1859.0 / 4104.0, f3);
        b5        = Add (b5,  -0.275,           f4);
        X    x5   = Add (*a_x, tau, b5);
        T    t5   = t + 0.5  * tau;
        auto f5   = a_rhs(x5,  t5);

        // Error estimate (NB: "f1" is not used here, this is correct):
        auto ee   = Mult(        1.0 / 360.0,   f0);
        ee        = Add (ee,  -128.0 / 4275.0,  f2);
        ee        = Add (ee, -2197.0 / 75240.0, f3);
        ee        = Add (ee,    0.02,           f4);
        ee        = Add (ee,     2.0 / 55.0,    f5);

        // NB: Max Abs Error (with "tau" factor), as a plain "double", since
        // the actual element type on which this max occurs, cannot be known
        // in advance:
        double maxAbsErr  = tauMag * MaxAbsRelError(ee, x5);
        assert(maxAbsErr >= 0.0);

        // XXX: If eps <= 0, do not use error / step control at all:
        if (a_eps <= 0.0 || maxAbsErr < a_eps)
        {
          // OK, the error bounds are satisfied, update "a_x". NB: "f1" is not
          // used, this is correct:
          auto b6 = Mult(       16.0 /   135.0, f0);
          b6      = Add (b6,  6656.0 / 12825.0, f2);
          b6      = Add (b6, 28561.0 / 56430.0, f3);
          b6      = Add(b6,    -0.18,           f4);
          b6      = Add(b6,      2.0 / 55.0,    f5);

          // This step is done:
          *a_x    = Add(*a_x, tau, b6);
          t       = tn;
          ok      = true;
          break;
        }
        else
        {
          // Error bound "eps" exceeded; we need to reduce the step and try
          // again. New step estimate:
          T tau1 = 0.9 * tau * std::pow(a_eps / maxAbsErr, 0.2);

          // However, we will not use "tau'" directly; better reduce the origi-
          // nal "tau" by powers of 2 until it becomes <= "tau1":
          while (std::fabs(tau) > std::fabs(tau1))
            tau *= 0.5;

          // And go for the next inner iteration with the new "tau"...
        }
      }
      // If all step reductions have been done, still to no avail:
      if (!ok)
        throw std::runtime_error
              ("StepRKF5: Required accuracy cannot be achieved");
    }
    // End of "t" loop: All Done!
  }
}
// End namespace SpaceBallistics
