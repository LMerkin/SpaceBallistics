// vim:ts=2:et
//===========================================================================//
//                            "PrecNutTest.cpp":                             //
//             Test for Long-Term Earth Precession and Nutations             //
//===========================================================================//
#include "SpaceBallistics/PhysForces/EarthRotationModel.hpp"
#include "SpaceBallistics/CoOrds/SpherPV.hpp"

int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;

  bool  nutAnglesTest = (argc == 2 && strcmp(argv[1], "-n") == 0);

  if (!nutAnglesTest)
  {
    //-----------------------------------------------------------------------//
    // Precession/Nutation Test over the full Precession Cycle from 2000.0:  //
    //-----------------------------------------------------------------------//
    // (The exact dates do not matter):
    constexpr Time_jyr From(2000.0);
    constexpr Time_jyr Step(1.0);
    constexpr Time_jyr To  = From + Time_jyr(26'000.0);

    for (Time_jyr t = From; t <= To; t += Step)
    {
      // ITRS orientation @ t:
      EarthRotationModel era{ t };

      // Get the PolarVec in GCRS:
      PosKV_GCRS<>  dynCelestPolarVec = era.GetGeoCDynEqFixZ();

      // Convert  the  "celestPolarVec" into Spherical Co-Ords:
      SpherPV<GCRS> dynCelestPolarSph(dynCelestPolarVec);

      // The Dynamic Equinox (of "t") vector extracted from the "era":    it's
      // Col0 of the PN matrix. Thus, it is the position of the Equinox of "t"
      // in the GCRS:
      PosKV_GCRS<>  dynEquinoxVec = era.GetGeoCDynEqFixX();
      SpherPV<GCRS> dynEquinoxSph(dynEquinoxVec);

      // Output:
      printf("%.1lf\t"                // Time (Year Number)
             "%.6lf\t"                // RA of the Dynamic Equinox
             "%.6lf\t%.6lf\t"         // "celestPolarVec" X  and Y
             "%.6lf\t%.6lf\n",        // "celestPolarSph" RA and Decl
             t.Magnitude(),
             dynEquinoxSph.GetAlpha()    .Magnitude(),
             dynCelestPolarVec.x()       .Magnitude(),
             dynCelestPolarVec.y()       .Magnitude(),
             dynCelestPolarSph.GetAlpha().Magnitude(),
             dynCelestPolarSph.GetDelta().Magnitude());
    }
  }
  else
  {
    //-----------------------------------------------------------------------//
    // Nutation Angles Test: Analytical vs DE440T:                           //
    //-----------------------------------------------------------------------//
    constexpr Time_jyr from = DE440T::Bits::FromY;
    constexpr Time_jyr to   = DE440T::Bits::ToY;
    constexpr Time_jyr step = 1.0_jyr / 12.0;   // 1m

    for (Time_jyr t = from; t < to; t += step)
    {
      auto nasAnalyt = EarthRotationModel::GetNutAnglesAnalyt(t);
      auto nasDE440T = EarthRotationModel::GetNutAnglesDE440T(t);

      Angle_arcSec dPhiAn = To_Angle_arcSec(nasAnalyt.first);
      Angle_arcSec dPhiDE = To_Angle_arcSec(nasDE440T.first);

      Angle_arcSec dEpsAn = To_Angle_arcSec(nasAnalyt.second);
      Angle_arcSec dEpsDE = To_Angle_arcSec(nasDE440T.second);

      printf("%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n",
             t.Magnitude(),
             dPhiAn.Magnitude(), dPhiDE.Magnitude(),
             dEpsAn.Magnitude(), dEpsDE.Magnitude());
    }
  }
  return 0;
}
