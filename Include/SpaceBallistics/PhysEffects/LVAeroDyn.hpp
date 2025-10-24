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
  // "Model" Enum Class:                                                     //
  //=========================================================================//
  enum class Model: int
  {
    M0 = 0,  // Bartholomeyev, Kopytov (from Titan II data)
    M1 = 1   // ChatGPT
  };

  //=========================================================================//
  // "cD": Drag Coeff:                                                       //
  //=========================================================================//
  inline double cD(double a_M, Angle a_AoA, Model a_model = Model::M0)
  {
    // FIXME: To avoid run-time exceptions, we set a VERY WIDE AoA limit here:
    // |AoA| <~ 23.0 deg; the correct one would probably be 5..8 deg. However,
    // large AoAs can only occur for upper stages, and for them the aerodynamic
    // effects are very small, anyway:
    //
    if (a_M < 0.0 || std::fabs(Sin(a_AoA)) > 0.4)
      throw std::invalid_argument
            ("cD: Invalid M=" + std::to_string(a_M) + " or AoA=" +
             std::to_string(double(a_AoA)));

    double cD0  = 0.0;

    switch (a_model)
    {
    case Model::M0:
      cD0 =
        (a_M <= 0.8)
        ? 0.29 :
        (a_M <= 1.068)
        ? a_M - 0.51
        : 0.091 + 0.5 / a_M;
      break;

    case Model::M1:
    {
      constexpr int    N = 12;
      constexpr double Ms [N]
      {0.0,  0.2,  0.4,  0.6,  0.8,  1.0,  1.2,  1.5,  2.0,  3.0,  5.0,  10.0};

      constexpr double cDs[N]
      {0.20, 0.20, 0.19, 0.18, 0.22, 0.60, 0.45, 0.40, 0.35, 0.35, 0.40, 0.45};

      // First, compute "cD0" for AoA=0. Direct linear search is optimal in
      // this case:
      for (int i = N-1; i >= 0; --i)
        if (Ms[i] <= a_M)
        {
          cD0 =
            (i == N-1)
            ? // There is no right value:
              cDs[i]
            :
              // Perform linear interpolation (XXX: the slope could be pre-
              // evaluated to improve efficiency):
              cDs[i] +
                ((cDs[i+1] - cDs[i]) / (Ms[i+1] - Ms[i])) * (a_M - Ms[i]);
          break;
        }
      break;
    }

    default: ;
    }

    // Now include the quadratic dependency on the AoA. The quadratic coeff is
    // typically 2..4; we use 3:
    assert(cD0 >= 0.0);
    return cD0 * (1.0 + 3.0 * Sqr(double(a_AoA)));
  }

  //=========================================================================//
  // "cL": Lift Coeff:                                                       //
  //=========================================================================//
  inline double cL(double a_M, Angle a_AoA, Model a_model = Model::M0)
  {
    // FIXME: Similar to "cD", we currently use an abnormally-large AoA limit:
    if (a_M < 0.0 || std::fabs(Sin(a_AoA)) > 0.4)
      throw std::invalid_argument
            ("cD: Invalid M=" + std::to_string(a_M) + " or AoA=" +
             std::to_string(double(a_AoA)));

    double cL1 = 0.0;

    switch (a_model)
    {
    case Model::M0:
      cL1 =
        (a_M <= 0.25)
        ? 2.8 :
        (a_M <=  1.1)
        ? 2.8  + 0.447 * (a_M - 0.25) :
        (a_M <=  1.6)
        ? 3.18 - 0.660 * (a_M - 1.1)  :
        (a_M <=  3.6)
        ? 2.85 + 0.350 * (a_M - 1.6)
        : 3.55;
      break;

    case Model::M1:
      if (a_M <= 0.8)
        cL1 = 2.0;
      else
      if (a_M <= 1.2)
      {
        double w  = 0.5 * (1.0 - TanH((a_M - 1.0) / 0.15));
        double d1 = SqRt(Max(1.0 - Sqr(a_M), 0.0) + 2.5e-4);
        double d2 = SqRt(Max(Sqr(a_M) - 1.0, 0.0) + 2.5e-4);
        cL1       = Min(3.0, 2.0 * (w / d1 + (1.0 - w) / d2));
      }
      else
        cL1       = 2.0 / SqRt(Sqr(a_M) - 1.0);
      break;

    default: ;
    }
    assert(cL1 >= 0.0);
    return cL1 * double(a_AoA);
  }
}
// End namespace SpaceBallistics::LVAeroDyn`
