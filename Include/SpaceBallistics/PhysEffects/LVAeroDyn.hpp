// vim:ts=2:et
//===========================================================================//
//               "SpaceBallistics/PhysEffects/LVAeroDyn.hpp":                //
//                                                                           //
//    Simplified Aprroximations for the AeroDynamic Drag and Lift Coeffs     //
//      (cD and cL, resp) as a function of the Mach Number (M) and the       //
//    Angle-of-Attack (AoA), for a Slender Cylindrical LV with a Pointed     //
//                                 Nose Cone                                 //
//===========================================================================//
#pragma once
#include <exception>
#include <cmath>

namespace SpaceBallistics::LVAeroDyn
{
  //=========================================================================//
  // "cD": Drag Coeff:                                                       //
  //=========================================================================//
  inline double cD(double a_M, Angle a_AoA)
  {
    constexpr int    N = 12;
    constexpr double Ms [N]
      {0.0,  0.2,  0.4,  0.6,  0.8,  1.0,  1.2,  1.5,  2.0,  3.0,  5.0,  10.0};

    constexpr double cDs[N]
      {0.20, 0.20, 0.19, 0.18, 0.22, 0.60, 0.45, 0.40, 0.35, 0.35, 0.40, 0.45};

    // First, compute "cD0" for AoA=0. Direct linear search is optimal in
    // this case:
    if (a_M < 0.0 || std::fabs(double(a_AoA) > 0.4)) // AoA Limit: +- 23.0 deg
      throw std::invalid_argument
            ("cD: Invalid M=" + std::to_string(a_M) + " or AoA=" +
             std::to_string(double(a_AoA)));

    double cD0 = 0.0;
    for (int i = N-1; i >= 0; --i)
      if (Ms[i] <= a_M)
      {
        cD0 =
          (i == N-1)
          ? // There is no right value:
            cDs[i]
          :
            // Perform linear interpolation (XXX: the slope could be pre-evalu-
            // ated):
            cDs[i] + ((cDs[i+1] - cDs[i]) / (Ms[i+1] - Ms[i])) * (a_M - Ms[i]);
        break;
      }
    assert(cD0 > 0.0);
    return cD0 * (1.0 + 3.0 * Sqr(double(a_AoA))); // Coeff: 2..4, we use 3
  }

  //=========================================================================//
  // "cL": Lift Coeff:                                                       //
  //=========================================================================//
  inline double cL(double a_M, Angle a_AoA)
  {
    // FIXME: TODO: XXX: Not implemented yet:
    return 0.0;
  }
}
// End namespace SpaceBallistics::LVAeroDyn`
