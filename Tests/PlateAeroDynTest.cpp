// vim:ts=2:et
//===========================================================================//
//                        "Tests/PlateAeroDynTest.cpp":                      //
//===========================================================================//
#include "SpaceBallistics/PhysEffects/PlateAeroDyn.hpp"
#include "SpaceBallistics/PhysEffects/EarthAtmosphereModel.h"
#include "SpaceBallistics/LVSC/LVSC.h"
#include "SpaceBallistics/CoOrds/EmbeddedCOS.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[])
{
  using namespace  SpaceBallistics;
  using namespace  std;
  namespace EAM  = EarthAtmosphereModel;
  using     ECOS = EmbeddedCOS<LVSC::Soyuz21b>;

  // The Angle of Attack and the normal to the Plate:
  if (argc != 2)
  {
    cerr << "PARAM: AngleOfAttack_deg" << endl;
    return 1;
  }
  Angle alpha = To_Angle(Angle_deg(atof(argv[1])));
  double sinAlpha = Sin(alpha);
  double cosAlpha = Cos(alpha);
  DimLessV<ECOS>        n { TT(), sinAlpha, cosAlpha, 0.0  };

  // Assume there is no wind:
  constexpr VelV<ECOS> w { TT(), Vel(0.0), Vel(0.0), Vel(0.0) };
  constexpr Area       S(1.0);

  // Plate velocity is towards (-X), in m/sec:
  //
  for (Vel Vx = Vel(10.0); Vx < Vel(1000.0); Vx += Vel(10.0))
  {
    VelV<ECOS>  v { TT(), -Vx,  Vel(0.0), Vel(0.0) };
    auto   [F, M, cR] =
      PlateAeroDyn(v, w, n, EAM::P0, EAM::Rho0, EAM::GammaAir, S);
    Force  R          =  F.EuclidNorm();
    double M2         =  Sqr(M);

    // Get the List and Drag coeffs from "cR":
    double cL         =  IsZero(R) ? 0.0 : cR * double(F.y() / R);
    double cD         =  IsZero(R) ? 0.0 : cR * double(F.x() / R);

    // And the approximate "cR" from a linearised model, for small "alpha",
    // in the supersonic case only:
    double cLappr     = (M > 1.4) ? double(alpha) * 4.0 / SqRt(M2 - 1.0) : cL;
    double cDappr     = (M > 1.4) ? double(alpha) * cLappr               : cD;

    cout << M      << '\t' << cR     << '\t'
         << cL     << '\t' << cD     << '\t'
         << cLappr << '\t' << cDappr << '\t'
         << F.x().Magnitude() / M2   << '\t'
         << F.y().Magnitude() / M2   << '\t'
         << F.z().Magnitude() / M2   << endl;
  }
  return 0;
}
