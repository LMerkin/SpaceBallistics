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

    // The Polar ITRS vector of unit length:
    constexpr PosKV<ITRS> TerrPolarVec { 0.0_km, 0.0_km, 1.0_km };

    for (Time_jyr t = From; t <= To; t += Step)
    {
      // ITRS orientation @ t:
      EarthRotationModel era{ t };

      // Convert the PolarVec to GCRS. For that, we have to convert "t" to TT,
      // although it does not matter here, since we only need  the Earth Axis
      // orientation:
      PosKV<GCRS> celestPolarVec = era.ToGCRS(TT{t}, TerrPolarVec);

      // Convert  the  "celestPolarVec" into Spherical Co-Ords:
      SpherPV<GCRS> celestPolarSph(celestPolarVec, VelKV<GCRS>());

      // The Dynamic Equinox vector extracted from the "era":  it's Col0 of the
      // PN matrix.
      // XXX: we don't use "GeoCDynEqFixCOS" here, so just use "Vector3D" with
      // DQ=LenK and COS=void, and construct "SpherPV" from it; ie it is a dir-
      // ectional vector of 1.0_km length:
      //
      Mtx33 const&   PN  = era.GetPN();
      Vector3D<LenK, void> equinox
        { 1.0_km * PN(0,0), 1.0_km * PN(1,0), 1.0_km * PN(2,0) };
      SpherPV<void> equinoxSph(equinox, VelKV<void>());

      // Output:
      printf("%.1lf\t"                // Time (Year Number)
             "%.6lf\t"                // RA of the Dynamic Equinox
             "%.6lf\t%.6lf\t"         // "celestPolarVec" X  and Y
             "%.6lf\t%.6lf\n",        // "celestPolarSph" RA and Decl
             t.Magnitude(),
             equinoxSph.GetAlpha()    .Magnitude(),
             celestPolarVec.x()       .Magnitude(),
             celestPolarVec.y()       .Magnitude(),
             celestPolarSph.GetAlpha().Magnitude(),
             celestPolarSph.GetDelta().Magnitude());
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
