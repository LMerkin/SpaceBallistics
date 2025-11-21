// vim:ts=2:et
//===========================================================================//
//                    "Src/Missions/Ascent2-Ret1.cpp":                       //
//                           Stage1 Return Params                            //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent2.h"

namespace SpaceBallistics
{
//===========================================================================//
// "SetStage1RetParams":                                                     //
//===========================================================================//
void Ascent2::SetStage1RetParams(LenK a_h, VelK a_V, Angle a_psi)
{
  // XXX: At the moment, only the RTLS mode is supported:
  assert(m_retMode1 == Stage1RetMode::RTLS);

  // Use the following pre-calibrated estimator for the PropMass (incl the Rem-
  // nants) required for the return and soft landing of Stage1.
  // XXX: Here the distance to the Landing Site is unknown, and should actually
  // be determined using the fixed-point method -- but for now, we just use a
  // constant:
  constexpr LenK x   = 75.0_km;
  double         h   = double(a_h   / 1.0_km);
  double         l   = double(x     / 1.0_km);
  double         V   = double(a_V   / VelK(1.0));
  double         psi = double(a_psi / 1.0_rad);

  m_propMassS1 =
    Mass
    (
      -183.450331400897539       + 1506.82402904732021*h
      -9.20379040331731879*h*h   + 1.04956247522381996*h*l
      +62.1640831399363663*V*h   - 7.53924156101525966*h*psi
      -1264.02219268784506*l     + 5.48619677522246185*l*l
      +112.539561322415437*V*l   - 2.42135702102538630*l*psi
      +2651.42176116741120*V     + 9712.42668440847956*V*V
      -42153.5700406520627*V*psi + 47037.4194317405781*psi
      -.359956980898263335e-1*psi*psi
    );
  // XXX: The above approximation is valid only for some range of the args, so
  // apply some reasonable constraints to the value computed:
  m_propMassS1 = std::max(30'000.0_kg, std::min(m_propMass1, 75'000.0_kg));

  if (Base::m_os != nullptr && Base::m_logLevel >= 3)
#   pragma omp critical(Output)
    *Base::m_os
      << "# Set PropMassS1=" << m_propMassS1.Magnitude()   << " kg: h="
      << a_h.Magnitude()     << " km, l=" << x.Magnitude() << " km, V="
      << a_V.Magnitude()     << " km/sec, psi="
      << To_Angle_deg(a_psi).Magnitude()  << " deg"        << std::endl;
}
}
// End namespace SpaceBallistics
