// vim:ts=2:et
//===========================================================================//
//                            "Tests/RKF5Test.cpp":                          //
//===========================================================================//
#include "SpaceBallistics/Maths/RKF5.hpp"
#include <tuple>
#include <iostream>

using namespace std;
using namespace SpaceBallistics;

namespace
{
  // ODE RHS:
  // x' = y
  // y' = x   ==> y'' =  x' = y
  // u' = v
  // v' = -u  ==> v'' = -u' = -v:
  //
  // Phase Space: s = [x, y, u, v]:
  using S = tuple<double, double, double, double>;

  S RHS(S const& a_s, double /*a_t*/)
  {
    return make_tuple
          (get<1>(a_s),  get<0>(a_s),
           get<3>(a_s), -get<2>(a_s));
  }

  // User Call-Back:
  bool CB(S const& a_x, double a_t)
  {
    // Expected solution:
    // x(t) = exp(t)
    // y(t) = exp(t)
    // u(t) = sin(t)
    // v(t) = cos(t)
    double e0 = fabs(get<0>(a_x) / exp(a_t) - 1.0);
    double e1 = fabs(get<1>(a_x) / exp(a_t) - 1.0);
    double e2 = fabs(get<2>(a_x) - sin(a_t));
    double e3 = fabs(get<3>(a_x) - cos(a_t));
    cout << a_t << "  " << e0 << "  " << e1 << "  " << e2 << "  " << e3 << endl;
    return true;
  }
}

int main()
{
  // Initial Conds:
  S s = make_tuple(1.0, 1.0, 0.0, 1.0);

  // Run the Solver for t=0..5, with the initial step "tau" and the required
  // accuracy "eps":
  constexpr double tau0 = 0.001;
  constexpr double eps  = 1e-12;

  RKF5(&s, 0.0, 5.0, RHS, tau0, eps, CB);
  return 0;
}
